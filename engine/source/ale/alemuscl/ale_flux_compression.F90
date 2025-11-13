!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      module ale_flux_compression_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
      subroutine ale_flux_compression(flux, ale_connect, x, ixq, flux_vois, n4_vois, itab, nv46, itrimat, segvar, &
                                      numelq, numnod, nqvois, nft, nel, nsegflu)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use i22bufbric_mod
      use i22tri_mod
      use alemuscl_mod , only:alemuscl_buffer
      use segvar_mod
      use ale_connectivity_mod
      use precision_mod , only : wp
      use element_mod , only : nixq
      use constant_mod , only : half, one, zero
      use ale_mod , only : ale
! -----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! -----------------------------------------------------------------------------------------------------------------------
          implicit none
! -----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! -----------------------------------------------------------------------------------------------------------------------
!#include "spmd_c.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer, intent(in) :: nv46
      integer,intent(in) :: numelq, numnod
      integer,intent(in) :: nsegflu !< number of boundary faces (inlet/outlet with EBCS)
      integer,intent(in) :: nqvois !< spmd parameter
      integer,intent(in) :: nft !< index shift (elem group)
      integer,intent(in) :: nel !< number of elems in the current group
      real(kind=wp), intent(out) :: flux(nv46, *)
      real(kind=wp), intent(in) :: x(3, numnod)
      integer, intent(in) :: ixq(nixq, numelq)
      real(kind=wp), intent(out) :: flux_vois(numelq+nqvois, nv46)
      integer, intent(out) :: n4_vois(numelq+nqvois,8)
      integer, intent(in) :: itab(numnod)
      integer, intent(in) :: itrimat
      type(t_segvar),intent(in) :: segvar
      type(t_ale_connectivity), intent(in) :: ale_connect
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: I, II, KK, JJ, IAD2, IAD3
      INTEGER :: NEIGHBOOR_LIST(NV46), FACE_NEIGHBOOR(NV46)
      real(kind=WP) :: ALPHAK
      real(kind=WP) :: YK, ZK
      real(kind=WP) :: YF, ZF
      real(kind=WP) :: kappa,norm_grad, corr, alpha_face,ny,nz,un_face,hard_coded_face
      INTEGER :: FACE_TO_NODE_LOCAL_ID(4, 2), NODEID1, NODEID2
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      kappa = ALE%INTER_CAPTURING%KAPPA
!!! Once for all, associate node local id to a face number
!!! Face 1
      FACE_TO_NODE_LOCAL_ID(1, 1) = 1 ; FACE_TO_NODE_LOCAL_ID(1, 2) = 2
!!! Face 2
      FACE_TO_NODE_LOCAL_ID(2, 1) = 2 ; FACE_TO_NODE_LOCAL_ID(2, 2) = 3
!!! Face 3
      FACE_TO_NODE_LOCAL_ID(3, 1) = 3 ; FACE_TO_NODE_LOCAL_ID(3, 2) = 4
!!! Face 4
      FACE_TO_NODE_LOCAL_ID(4, 1) = 4 ; FACE_TO_NODE_LOCAL_ID(4, 2) = 1
!!! First of all, compute gradient for alpha
      DO I = 1,NEL
        II = I + NFT
        IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
        !!!centroid element
        YK = ALEMUSCL_Buffer%ELCENTER(II,2) ;
        ZK = ALEMUSCL_Buffer%ELCENTER(II,3)
        !!! Neighbors
        DO KK = 1, NV46
          !!! Only for outgoing fluxes
          IF (FLUX(KK, II) > ZERO) THEN
            !!! Storing neighbor indexes
            NEIGHBOOR_LIST(KK) = ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1)
            FACE_NEIGHBOOR(KK) = KK
            IF (NEIGHBOOR_LIST(KK) <= 0) THEN
              IF(NEIGHBOOR_LIST(KK)==0)NEIGHBOOR_LIST(KK) = II
              !case <0 is for eBCS. -NEIGHBOR_LIST is then the segment number
            ELSEIF (NEIGHBOOR_LIST(KK) <= NUMELQ) THEN
              IAD3 = ALE_CONNECT%ee_connect%iad_connect(NEIGHBOOR_LIST(KK))
              !!! Store the face number to which II and NEIGHBOR_LIST(KK) are adjacent
              DO JJ = 1, NV46
                IF (ALE_CONNECT%ee_connect%connected(IAD3 + JJ - 1) == II) THEN
                  FACE_NEIGHBOOR(KK) = JJ
                ENDIF
              ENDDO  ! JJ = 1, NV46
            ENDIF

            NODEID1 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 1), II)
            NODEID2 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 2), II)

            YF = HALF * (X(2, NODEID1) + X(2, NODEID2))
            ZF = HALF * (X(3, NODEID1) + X(3, NODEID2))

            !!! Reconstruct second order value for ALPHA(II) on the face
            ALPHAK = ALEMUSCL_Buffer%VOLUME_FRACTION(II,ITRIMAT) &
                 + ALEMUSCL_Buffer%GRAD(II,2,ITRIMAT) * (YF - YK)&
                 + ALEMUSCL_Buffer%GRAD(II,3,ITRIMAT) * (ZF - ZK)

            !!! Partial volume flux is then computed as:
            FLUX(KK, II) = ALPHAK * FLUX(KK, II)





            !      --- dans la boucle sur les faces (KK,II) ---
            !      calcul du |∇α| et du α moyen sur la face
                   ny = HALF*( ALEMUSCL_Buffer%GRAD(II,2,ITRIMAT)+ ALEMUSCL_Buffer%GRAD(NEIGHBOOR_LIST(KK),2,ITRIMAT) )
                   nz = HALF*( ALEMUSCL_Buffer%GRAD(II,3,ITRIMAT)+ ALEMUSCL_Buffer%GRAD(NEIGHBOOR_LIST(KK),3,ITRIMAT) )
                   norm_grad = sqrt(ny*ny + nz*nz)

                   alpha_face = HALF*( ALEMUSCL_Buffer%VOLUME_FRACTION(II,ITRIMAT) &
                              + ALEMUSCL_Buffer%VOLUME_FRACTION(NEIGHBOOR_LIST(KK),ITRIMAT) )

            !      --- on applique la correction uniquement sur les "rides" ---
                   IF (norm_grad > ZERO .and. alpha_face > 0.01 .and. alpha_face < 0.99) THEN
            !         vitesse normale à la face (déjà stockée dans FLUX global)
            !          hard_coded_face = 0.0025  !not stored in ALEMUSCL_Buffer
            !          un_face = -(FLUX(KK,II)) / hard_coded_face   ! m/s
                      corr = -KAPPA * (-FLUX(KK,II)) * alpha_face*(1.0 - alpha_face)          ! m³/s
                      corr = max( -0.05*abs(FLUX(KK,II)), min( 0.05*abs(FLUX(KK,II)), corr ) )
            !         on retire ce flux à la face "sortante" et on l’ajoute à la face "entrante"
                      FLUX(KK,II) = FLUX(KK,II) + corr
                       IF (NEIGHBOOR_LIST(KK) > 0 .and. NEIGHBOOR_LIST(KK) <= NUMELQ) THEN
                          FLUX(FACE_NEIGHBOOR(KK), NEIGHBOOR_LIST(KK)) = &
                             FLUX(FACE_NEIGHBOOR(KK), NEIGHBOOR_LIST(KK)) - corr
                       ENDIF
                   ENDIF




            IF (NEIGHBOOR_LIST(KK) > 0)THEN
              IF (NEIGHBOOR_LIST(KK) <= NUMELQ) THEN
                !!! The opposite of the flux goes to the neighbor
                FLUX(FACE_NEIGHBOOR(KK), NEIGHBOOR_LIST(KK)) = -FLUX(KK, II)
              ELSE
                !!! cf. ALE51_ANTIDIFF3
                FLUX_VOIS(II, KK) = FLUX(KK, II)
                N4_VOIS(II, 1) = ITAB(IXQ(2, II))
                N4_VOIS(II, 2) = ITAB(IXQ(3, II))
                N4_VOIS(II, 3) = ITAB(IXQ(4, II))
                N4_VOIS(II, 4) = ITAB(IXQ(5, II))
              ENDIF
            ENDIF
          ENDIF  ! (FLUX(KK, II) > ZERO)
        ENDDO  ! KK = 1, NV46
      ENDDO  ! I = 1,NEL

      !-----------------------------------------------
      !         incoming volume flux (EBCS)
      !-----------------------------------------------
      IF(NSEGFLU > 0)THEN
        DO I = 1,NEL
          II = I + NFT
          IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
          DO KK=1,4
            IF(FLUX(KK,II) < ZERO .AND. ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1) < 0)THEN
              FLUX(KK,II) = SEGVAR%PHASE_ALPHA(ITRIMAT,-ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1))*FLUX(KK,II)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------
      END SUBROUTINE ale_flux_compression
      END MODULE ale_flux_compression_mod
