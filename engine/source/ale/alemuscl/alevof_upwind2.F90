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
!||====================================================================
!||    alevof_upwind2               ../engine/source/ale/alemuscl/alevof_upwind2.F90
!||--- called by ------------------------------------------------------
!||    afluxt                       ../engine/source/ale/ale51/afluxt.F
!||--- uses       -----------------------------------------------------
!||    ale_connectivity_mod         ../common_source/modules/ale/ale_connectivity_mod.F
!||    ale_mod                      ../common_source/modules/ale/ale_mod.F
!||    constant_mod                 ../common_source/modules/constant_mod.F
!||    element_mod                  ../common_source/modules/elements/element_mod.F90
!||    precision_mod                ../common_source/modules/precision_mod.F90
!||    segvar_mod                   ../engine/share/modules/segvar_mod.F
!||====================================================================
!! \brief VOF-PLIC upwind for 2D quads — 2 phases only
!! \details Computes sub-material volume fluxes using the PLIC reconstruction:
!!          - For PURE cells (alpha=0 or 1): flux_mat(:,:,1) = alpha * flux_total,
!!            flux_mat(:,:,2) = (1-alpha) * flux_total  (already set by caller)
!!          - For MIXED cells: the wet fraction Swet(KK)/edge_length determines
!!            the phase-1 share of the outgoing flux through each face KK.
!!            Phase 2 gets the complement.
!!          - For INCOMING fluxes: the upwind cell's alpha determines the phase split.
!!          - Phases 3..trimat get zero flux.
       MODULE ALEVOF_UPWIND2_MOD
         IMPLICIT NONE
       CONTAINS
       SUBROUTINE ALEVOF_UPWIND2(flux_mat,flux_save, ALE_CONNECT, X, IXQ, flux_vois_mat, &
           NV46, trimat, SEGVAR,s_flux,s_flux_vois, &
           NUMELQ,NUMNOD,NSEGFLU,NEL,NFT)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE SEGVAR_MOD
      USE ALE_CONNECTIVITY_MOD
      use element_mod , only : nixq
      use precision_mod , only : WP
      use constant_mod , only : HALF, ZERO, ONE, EM20
      use ale_mod , only : ale
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Dummy Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer, intent(in) :: s_flux       !< first dimension of FLUX array
      integer, intent(in) :: s_flux_vois  !< first dimension of flux_vois_mat
      integer, intent(in) :: trimat       !< number of materials (must be >= 2)
      INTEGER, INTENT(IN) :: NV46
      real(kind=WP), dimension(s_flux,nv46), intent(in) :: flux_save
      real(kind=WP), dimension(s_flux,NV46,trimat), INTENT(INOUT) :: flux_mat
      real(kind=WP), INTENT(IN) :: X(3, NUMNOD)
      INTEGER, INTENT(IN) :: IXQ(NIXQ, NUMELQ)
      real(kind=WP), INTENT(INOUT) :: flux_vois_mat(s_flux_vois, NV46,trimat)
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
      INTEGER :: I, II, KK, JJ, IAD2
      INTEGER :: NEIGH, FACE_NEIGH
      INTEGER :: ICELL
      INTEGER :: FACE_TO_NODE_LOCAL_ID(4, 2), NODEID1, NODEID2
      real(kind=WP) :: ALPHA_I, ALPHA_FACE
      real(kind=WP) :: EDGE_LEN, WET_FRAC
      real(kind=WP) :: FLUX_TOTAL
      real(kind=WP) :: Y1, Z1, Y2, Z2
      real(kind=WP), parameter :: TOL_PURE = 1.0e-3_WP  ! must match tol1 in ale51_vof_reconstruction2
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

      ! Face-to-node connectivity for QUAD4
      FACE_TO_NODE_LOCAL_ID(1, 1) = 1 ; FACE_TO_NODE_LOCAL_ID(1, 2) = 2
      FACE_TO_NODE_LOCAL_ID(2, 1) = 2 ; FACE_TO_NODE_LOCAL_ID(2, 2) = 3
      FACE_TO_NODE_LOCAL_ID(3, 1) = 3 ; FACE_TO_NODE_LOCAL_ID(3, 2) = 4
      FACE_TO_NODE_LOCAL_ID(4, 1) = 4 ; FACE_TO_NODE_LOCAL_ID(4, 2) = 1

      ! -----------------------------------------------
      ! NOTE: The default initialization (phases 1/2 = alpha-weighted, phases 3..trimat = 0)
      ! is done by the caller (afluxt.F) BEFORE the group loop.
      ! It CANNOT be done here because this routine is called per-group,
      ! and the init pass would overwrite neighbor flux values written by previous groups.
      ! -----------------------------------------------

      ! -----------------------------------------------
      ! OUTGOING FLUXES (flux_save > 0) : VOF phase split for phases 1 and 2
      ! -----------------------------------------------
      DO I = 1, NEL
        II = I + NFT
        IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
        ICELL = ALE%VOF%cell_data%mixed_cell_id(II)  ! 0 if pure, >0 if mixed

        DO KK = 1, NV46
          FLUX_TOTAL = flux_save(II, KK)

          IF(FLUX_TOTAL > ZERO) THEN

            IF( ICELL > 0) THEN
              ! --- MIXED cell: use Swet for phase split ---
              NODEID1 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 1), II)
              NODEID2 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 2), II)
              Y1 = X(2, NODEID1) ; Z1 = X(3, NODEID1)
              Y2 = X(2, NODEID2) ; Z2 = X(3, NODEID2)
              EDGE_LEN = SQRT((Y2 - Y1)**2 + (Z2 - Z1)**2)
              WET_FRAC = ALE%VOF%cell_data%Swet(KK, II) / (EDGE_LEN + EM20)
              WET_FRAC = MAX(ZERO, MIN(ONE, WET_FRAC))
              ! Phase 1 flux
              flux_mat(II, KK, 1) = WET_FRAC * FLUX_TOTAL
              ! Phase 2 flux (complement)
              flux_mat(II, KK, 2) = (ONE - WET_FRAC) * FLUX_TOTAL
            ELSE
              ! --- PURE cell: use cell alpha for phase split ---
              ALPHA_I = ALE%VOF%cell_data%ALPHA(II)
              ! Clip near-zero / near-one alpha to prevent micro-flux leakage
              ! that would spuriously create mixed cells in neighbors.
              IF(ALPHA_I < TOL_PURE) THEN
                ALPHA_I = ZERO
              ELSEIF(ALPHA_I > ONE - TOL_PURE) THEN
                ALPHA_I = ONE
              ENDIF
              flux_mat(II, KK, 1) = ALPHA_I * FLUX_TOTAL
              flux_mat(II, KK, 2) = (ONE - ALPHA_I) * FLUX_TOTAL
            ENDIF
            ! Zero out phases 3..trimat
            DO JJ = 3, trimat
              flux_mat(II, KK, JJ) = ZERO
            ENDDO

            ! --- Apply to neighbor (incoming flux = -outgoing flux) ---
            NEIGH = ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1)

            IF(NEIGH > 0) THEN
              IF(NEIGH <= NUMELQ) THEN
                ! Local neighbor: find the face number on the neighbor side
                FACE_NEIGH = ALE_CONNECT%ee_connect%IFACE2(IAD2 + KK - 1)
                IF(FACE_NEIGH > 0) THEN
                  ! The opposite flux goes to the neighbor
                  flux_mat(NEIGH, FACE_NEIGH, 1) = -flux_mat(II, KK, 1)
                  flux_mat(NEIGH, FACE_NEIGH, 2) = -flux_mat(II, KK, 2)
                  DO JJ = 3, trimat
                    flux_mat(NEIGH, FACE_NEIGH, JJ) = ZERO
                  ENDDO
                ENDIF
              ELSE
                ! Remote neighbor (SPMD): store in vois buffer
                flux_vois_mat(II, KK, 1) = flux_mat(II, KK, 1)
                flux_vois_mat(II, KK, 2) = flux_mat(II, KK, 2)
                DO JJ = 3, trimat
                  flux_vois_mat(II, KK, JJ) = ZERO
                ENDDO
              ENDIF
            ENDIF
            ! NEIGH==0 : free boundary, nothing to do
            ! NEIGH<0  : eBCS segment, handled in incoming flux section below

          ENDIF  ! FLUX_TOTAL > ZERO

        ENDDO  ! KK
      ENDDO  ! I

      ! -----------------------------------------------
      ! INCOMING FLUXES BY EBCS (boundary conditions)
      ! -----------------------------------------------
      IF(NSEGFLU > 0) THEN
        DO I = 1, NEL
          II = I + NFT
          IAD2 = ALE_CONNECT%ee_connect%iad_connect(II)
          DO KK = 1, 4
            IF(flux_save(II, KK) < ZERO .AND. &
               ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1) < 0) THEN
              FLUX_TOTAL = flux_save(II, KK)
              ! Use boundary phase fraction from SEGVAR
              ALPHA_FACE = SEGVAR%PHASE_ALPHA(1, -ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1))
              flux_mat(II, KK, 1) = ALPHA_FACE * FLUX_TOTAL
              flux_mat(II, KK, 2) = (ONE - ALPHA_FACE) * FLUX_TOTAL
              DO JJ = 3, trimat
                flux_mat(II, KK, JJ) = ZERO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF


      END SUBROUTINE ALEVOF_UPWIND2
      END MODULE ALEVOF_UPWIND2_MOD
