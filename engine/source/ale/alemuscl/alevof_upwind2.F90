!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
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
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief VOF method for 2D quads
!! \details      1 - build volume fraction at cell edge using the VOF method.
!! \details      2 - upwind this value on the edge and store it in the volumetric flux
!! \details if needed, more details can be added here
       MODULE ALEVOF_UPWIND2_MOD
         IMPLICIT NONE
       CONTAINS
       SUBROUTINE ALEVOF_UPWIND2(flux_mat,flux_save, ALE_CONNECT, X, IXQ, flux_vois_mat, &
           NV46, trimat, SEGVAR,s_flux,s_flux_vois, &
           NUMELQ,NUMNOD,NSEGFLU,NEL,NFT)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE I22BUFBRIC_MOD
      USE I22TRI_MOD
      USE ALEMUSCL_MOD , only:ALEMUSCL_Buffer
      USE SEGVAR_MOD
      USE ALE_CONNECTIVITY_MOD
      use element_mod , only :nixq
      use precision_mod , only : WP
      use constant_mod , only : half,  zero
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
      implicit none 
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Include Files
! ----------------------------------------------------------------------------------------------------------------------
#include "spmd.inc"
!#include "vect01_c.inc"
!#include "com04_c.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      integer, intent(in) :: s_flux !< first dimension of FLUX array
      integer, intent(in) :: s_flux_vois !< first dimension of flux_vois_mat
      integer, intent(in) :: trimat !< number of material      
      INTEGER, INTENT(IN) :: NV46
      real(kind=WP), dimension(s_flux,nv46), intent(in) :: flux_save
      real(kind=WP), dimension(s_flux,NV46,trimat), INTENT(OUT) :: flux_mat
      real(kind=WP), INTENT(IN) :: X(3, NUMNOD)
      INTEGER, INTENT(IN) :: IXQ(NIXQ, NUMELQ)
      real(kind=WP), INTENT(OUT) :: flux_vois_mat(s_flux_vois, NV46,trimat)
      TYPE(t_segvar),INTENT(IN) :: SEGVAR
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
      INTEGER,INTENT(IN) :: NUMELQ
      INTEGER,INTENT(IN) :: NUMNOD
      INTEGER,INTENT(IN) :: NSEGFLU
      INTEGER,INTENT(IN) :: NEL
      INTEGER,INTENT(IN) :: NFT
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: I, II, KK, JJ, IAD2, IAD3
      integer :: itrimat
      INTEGER :: NEIGHBOOR_LIST(NV46), FACE_NEIGHBOOR(NV46)
      real(kind=WP) :: ALPHAK
      real(kind=WP) :: YK, ZK
      real(kind=WP) :: YF, ZF
      INTEGER :: FACE_TO_NODE_LOCAL_ID(4, 2), NODEID1, NODEID2
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
!!! Once for all, associate node local id to a face number
      FACE_TO_NODE_LOCAL_ID(1, 1) = 1 ; FACE_TO_NODE_LOCAL_ID(1, 2) = 2  !!! Face 1
      FACE_TO_NODE_LOCAL_ID(2, 1) = 2 ; FACE_TO_NODE_LOCAL_ID(2, 2) = 3  !!! Face 2
      FACE_TO_NODE_LOCAL_ID(3, 1) = 3 ; FACE_TO_NODE_LOCAL_ID(3, 2) = 4  !!! Face 3
      FACE_TO_NODE_LOCAL_ID(4, 1) = 4 ; FACE_TO_NODE_LOCAL_ID(4, 2) = 1  !!! Face 4
!!! First of all, compute gradient for alpha
      DO I = 1, NEL
        II = I + NFT
        IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
        !!!centroid element
        YK = ALEMUSCL_Buffer%ELCENTER(II,2) ; 
        ZK = ALEMUSCL_Buffer%ELCENTER(II,3)
        !!! Neighbors
        DO KK = 1, NV46
          !!! Only for outgoing fluxes
          IF (flux_save(II,KK) > ZERO) THEN
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
            do itrimat=1,trimat
              ALPHAK = ALEMUSCL_Buffer%VOLUME_FRACTION(II,ITRIMAT) &
                      + ALEMUSCL_Buffer%GRAD(II,2,ITRIMAT) * (YF - YK) &
                      + ALEMUSCL_Buffer%GRAD(II,3,ITRIMAT) * (ZF - ZK)

              !!! Partial volume flux is then computed as:
              flux_mat(II,KK,itrimat) = ALPHAK * flux_save(II,KK)
            enddo

            IF (NEIGHBOOR_LIST(KK) > 0)THEN
              IF (NEIGHBOOR_LIST(KK) <= NUMELQ) THEN
                do itrimat=1,trimat                
                  !!! The opposite of the flux goes to the neighbor
                  flux_mat(NEIGHBOOR_LIST(KK),FACE_NEIGHBOOR(KK),itrimat) = -flux_mat(II,KK,itrimat)
                enddo
              ELSE
                do itrimat=1,trimat
                  flux_vois_mat(II, KK,itrimat) = flux_mat(II,KK,itrimat)
                enddo
              ENDIF
            ENDIF
          ENDIF  ! (FLUX(II,KK) > ZERO)
        ENDDO  ! KK = 1, NV46
      ENDDO  ! I = 1, NEL

!-----------------------------------------------
!         incoming flux by EBCS
!-----------------------------------------------
      IF(NSEGFLU > 0)THEN
        DO I = 1, NEL
          II = I + NFT
          IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
          DO KK=1,4            
            IF(flux_save(II,KK) < ZERO .AND. ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1) < 0)THEN
              do itrimat=1,trimat
                flux_mat(II,KK,itrimat) = SEGVAR%PHASE_ALPHA(ITRIMAT,-ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1))* &
                                         flux_save(II,KK)
              enddo
            ENDIF            
          ENDDO
        ENDDO
      ENDIF


      END SUBROUTINE ALEVOF_UPWIND2
      END MODULE ALEVOF_UPWIND2_MOD
