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
           NUMELQ,NUMNOD,NSEGFLU,NEL,NFT,DT1)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE SEGVAR_MOD
      USE ALE_CONNECTIVITY_MOD
      use element_mod , only : nixq
      use precision_mod , only : WP
      use constant_mod , only : HALF, ZERO, ONE, EM20
      use ale_mod , only : ale
      use polygon_mod , only : polygon_, polygon_list_, polygon_create, polygon_addpoint, polygon_point_
      use polygon_clipping_mod , only : Clipping_Weiler_Atherton
      use ale51_vof_reconstruction2_mod , only : CLIPPED_AREA_QUAD4
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
      real(kind=WP), INTENT(IN) :: DT1 !< time step
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
      INTEGER :: I, II, KK, JJ, IAD2, IERR, ISTATUS
      INTEGER :: NEIGH, FACE_NEIGH
      INTEGER :: ICELL
      INTEGER :: FACE_TO_NODE_LOCAL_ID(4, 2), NODEID1, NODEID2
      INTEGER :: ITRIMAT, IELEM
      real(kind=WP) :: ALPHA_I, ALPHA_FACE
      real(kind=WP) :: EDGE_LEN, INV_EDGE_LEN, WET_FRAC
      real(kind=WP) :: FLUX_TOTAL
      real(kind=WP) :: Y1, Z1, Y2, Z2
      real(kind=WP) :: NY, NZ, DY, DZ  !< inward normal components and displacement
      real(kind=WP) :: D_PLIC         !< PLIC plane distance
      real(kind=WP) :: NN_PLIC(3)     !< PLIC plane normal
      real(kind=WP) :: PTS_SWEPT(3,4) !< swept polygon vertices in PTS format
      real(kind=WP) :: CLIPPED_SWEPT  !< area of swept polygon clipped by PLIC half-plane
      real(kind=WP), parameter :: TOL_PURE = 1.0e-3_WP  ! must match tol1 in ale51_vof_reconstruction2
      real(kind=WP) :: ALPHA_MAIN !< main phase fraction in swept region
      real(kind=WP) :: ALPHA_REMAIN !< remaining ratio (1 - alpha_main)
      real(kind=WP) :: SUM_VF !< sum of non-main phase fractions
      real(kind=WP) :: VFRAC(4)
      type(polygon_) :: swept_polygon
      type(polygon_point_) :: pt
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      CALL polygon_create(swept_polygon, 4) ! swept polygon (advected volume) has always 4 nodes

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
            !outgoing flux

            IF( ICELL > 0) THEN
              ! --- MIXED cell: build swept polygon ---

              ! CONSTRUCT SWEPT POLYGON
              ! The quad is CCW in (y,z) plane.
              ! Inward normal of edge [N1,N2] is: n_in = (-(Z2-Z1), +(Y2-Y1)) / |edge|
              ! Swept depth = FLUX_TOTAL / EDGE_LEN (volume flux = area in 2D = depth * edge_len)
              NODEID1 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 1), II)
              NODEID2 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(KK, 2), II)
              Y1 = X(2, NODEID1) ; Z1 = X(3, NODEID1)
              Y2 = X(2, NODEID2) ; Z2 = X(3, NODEID2)
              EDGE_LEN = SQRT((Y2 - Y1)**2 + (Z2 - Z1)**2)
              INV_EDGE_LEN = ONE / MAX(EDGE_LEN, EM20)
              ! Inward unit normal (CCW connectivity → rotate edge vector -90°)
              NY = -(Z2 - Z1) * INV_EDGE_LEN
              NZ =  (Y2 - Y1) * INV_EDGE_LEN
              ! Inward displacement = (flux_rate * DT1) / edge_length = swept depth
              DY = NY * FLUX_TOTAL * DT1 * INV_EDGE_LEN
              DZ = NZ * FLUX_TOTAL * DT1 * INV_EDGE_LEN
              ! Build swept quad: N1 → N2 → N2+d → N1+d (CCW)
              swept_polygon%numpoint = 0
              pt%y = Y1 ;      pt%z = Z1
              IERR = polygon_addpoint(swept_polygon, pt)
              pt%y = Y2 ;      pt%z = Z2
              IERR = polygon_addpoint(swept_polygon, pt)
              pt%y = Y2 + DY ; pt%z = Z2 + DZ
              IERR = polygon_addpoint(swept_polygon, pt)
              pt%y = Y1 + DY ; pt%z = Z1 + DZ
              IERR = polygon_addpoint(swept_polygon, pt)
              swept_polygon%area = FLUX_TOTAL*DT1
              ! SWEPT POLYGON CONSTRUCTED

              !CLIPPING : SWEPT POLYGON ∩ PLIC half-plane {x | n.x >= d}
              ! Clip the swept quad (4 vertices) by the PLIC plane using Sutherland-Hodgman
              D_PLIC = ALE%VOF%cell_data%d(ICELL)
              NN_PLIC(2) = ALE%VOF%cell_data%n(2, ICELL)
              NN_PLIC(3) = ALE%VOF%cell_data%n(3, ICELL)
              ! Build PTS_SWEPT(3,4) from swept_polygon points
              PTS_SWEPT(1, 1:4) = ZERO
              DO JJ = 1, 4
                PTS_SWEPT(2, JJ) = swept_polygon%point(JJ)%y
                PTS_SWEPT(3, JJ) = swept_polygon%point(JJ)%z
              ENDDO
              CALL CLIPPED_AREA_QUAD4(PTS_SWEPT, NN_PLIC, D_PLIC, CLIPPED_SWEPT)

              alpha_main = CLIPPED_SWEPT / MAX(swept_polygon%area, EM20)
              alpha_main = MIN(alpha_main, ONE)  ! cap at 1 (geometric safety)
              ITRIMAT =  ALE%VOF%cell_data%IPHASE(II)

              VFRAC(1:4) = ALE%VOF%cell_data%ALPHA(II,1:4)  ! phase fractions in the cell (sum to 1)
              VFRAC(ITRIMAT) = ZERO
              SUM_VF = SUM(VFRAC(1:4))
              IF(SUM_VF > EM20) THEN
                VFRAC(1:4) = VFRAC(1:4) / SUM_VF
              ELSE
                VFRAC(1:4) = ZERO
              ENDIF
              ALPHA_REMAIN = ONE - alpha_main
              flux_mat(II, KK, 1) = VFRAC(1)*ALPHA_REMAIN * FLUX_TOTAL
              flux_mat(II, KK, 2) = VFRAC(2)*ALPHA_REMAIN * FLUX_TOTAL
              flux_mat(II, KK, 3) = VFRAC(3)*ALPHA_REMAIN * FLUX_TOTAL
              flux_mat(II, KK, 4) = VFRAC(4)*ALPHA_REMAIN * FLUX_TOTAL

              flux_mat(II, KK, ITRIMAT) = alpha_main * FLUX_TOTAL

            ELSE
              ! --- PURE cell => 1st order UPWIND ---
              flux_mat(II, KK, 1) = ALE%VOF%cell_data%ALPHA(II,1) * FLUX_TOTAL
              flux_mat(II, KK, 2) = ALE%VOF%cell_data%ALPHA(II,2) * FLUX_TOTAL
              flux_mat(II, KK, 3) = ALE%VOF%cell_data%ALPHA(II,3) * FLUX_TOTAL
              flux_mat(II, KK, 4) = ALE%VOF%cell_data%ALPHA(II,4) * FLUX_TOTAL
            ENDIF

            ! --- Apply to neighbor (incoming flux = -outgoing flux) ---
            NEIGH = ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1)

            IF(NEIGH > 0) THEN
              IF(NEIGH <= NUMELQ) THEN
                ! Local neighbor: find the face number on the neighbor side
                FACE_NEIGH = ALE_CONNECT%ee_connect%IFACE2(IAD2 + KK - 1)
                IF(FACE_NEIGH > 0) THEN
                  ! The opposite flux goes to the neighbor
                  DO JJ = 1, trimat
                    flux_mat(NEIGH, FACE_NEIGH, JJ) = -flux_mat(II, KK, JJ)
                  ENDDO
                ENDIF
              ELSE
                ! Remote neighbor (SPMD): store in vois buffer
                DO JJ = 1, trimat
                  flux_vois_mat(II, KK, JJ) = flux_mat(II, KK, JJ)
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
              IELEM = -ALE_CONNECT%ee_connect%connected(IAD2 + KK - 1)
              flux_mat(II, KK, 1) = SEGVAR%PHASE_ALPHA(1, IELEM) * FLUX_TOTAL
              flux_mat(II, KK, 2) = SEGVAR%PHASE_ALPHA(2, IELEM) * FLUX_TOTAL
              flux_mat(II, KK, 3) = SEGVAR%PHASE_ALPHA(3, IELEM) * FLUX_TOTAL
              flux_mat(II, KK, 4) = SEGVAR%PHASE_ALPHA(4, IELEM) * FLUX_TOTAL
            ENDIF
          ENDDO
        ENDDO
      ENDIF


      END SUBROUTINE ALEVOF_UPWIND2
      END MODULE ALEVOF_UPWIND2_MOD
