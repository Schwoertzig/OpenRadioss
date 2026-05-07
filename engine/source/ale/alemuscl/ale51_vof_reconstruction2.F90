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
      MODULE ALE51_VOF_RECONSTRUCTION2_MOD
        IMPLICIT NONE
      CONTAINS
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief This cell is Volume of FLuid Methods for 2D quads and 2 phases and unstructured meshes
!! \details step 1: compute volume fraction at cell center and store it in ALE%VOF%cell_data%ALPHA(IE)
!! \details step 2: compute volume fraction at cell edge and store it in ALE%VOF%cell_data%ALPHA_F(JJ, IE) where JJ is the local edge
!! \details step 3: compute normalized gradient of volume fraction for each cell(interface normal) and store it in ALE%VOF%cell_data%n(1:3, IE) where n(1) is always 0 in 2D
!! \details step 4: get d value (n-projection) for each point composing the quad, retain dmin and dmax. Find d value such as intersection is alpha


      SUBROUTINE ALE51_VOF_RECONSTRUCTION2(IPARG   ,ELBUF_TAB      ,ALE_CONNECT, NIXQ,IXQ, NUMNOD,X, &
                                           NERCVOIS ,NESDVOIS,LERCVOIS,LESDVOIS   ,LENCOM, ITASK, &
                                           SEGVAR  ,timers , &
                                           NPARG,NGROUP,NUMELQ,NUMELS,TRIMAT,NQVOIS,NSVOIS,NSPMD)

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Comments
! ----------------------------------------------------------------------------------------------------------------------
!
!    For each mixed cell :
!       1. Get 4 quad node coordinates PTS(2:3, 1:4) and interface normal NN(2:3)
!       2. Compute projections D(1:4) = NN . PTS for each node
!       3. Dmin = min(D), Dmax = max(D)  → bounds for Brent search
!       4. A_TARGET = alpha * AREA_CELL
!       5. Find d in [Dmin, Dmax] such that CLIPPED_AREA_QUAD4(PTS, NN, d) = A_TARGET
!          using Brent's method (guaranteed convergence, ~10 iterations)
!       6. Store d → ALE%VOF%cell_data%d(ICELL)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE INITBUF_MOD
      USE ELBUFDEF_MOD
      USE TRIMAT_MOD
      USE SEGVAR_MOD
      USE ALE_CONNECTIVITY_MOD
      USE MULTIMAT_PARAM_MOD , ONLY : M51_N0PHAS, M51_NVPHAS
      use timer_mod
      use precision_mod , only : WP      
      use spmd_exch_min_max_mod , only : spmd_exch_min_max
      use constant_mod , only : ZERO, ONE, EP20, EM20, FOURTH, HALF, EM12
      use ale_mod , only : ale
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit None
! ----------------------------------------------------------------------------------------------------------------------
      IMPLICIT     NONE
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
#include "task_c.inc"
#include "spmd.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Dummy Arguments
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: ITASK
      INTEGER IPARG(NPARG,NGROUP)
      TYPE(ELBUF_STRUCT_), TARGET, DIMENSION(NGROUP) :: ELBUF_TAB
      INTEGER :: LENCOM, NERCVOIS(*),NESDVOIS(*),LERCVOIS(*),LESDVOIS(*)
      TYPE(t_segvar) :: SEGVAR
      TYPE(t_ale_connectivity), INTENT(IN) :: ALE_CONNECT
      type(timer_), intent(inout) :: timers
      INTEGER,INTENT(IN) :: NPARG, NGROUP, NUMELQ, NUMELS,TRIMAT,NQVOIS,NSVOIS,NSPMD
      integer,intent(in) :: nixq
      INTEGER,intent(in) ::  IXQ(NIXQ, NUMELQ)
      INTEGER,INTENT(IN) :: NUMNOD
      real(kind=WP), intent(in) :: X(3, NUMNOD)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: NG
      INTEGER :: ITRIMAT
      real(kind=WP), DIMENSION(:), POINTER :: VOLG, VOLP, UVAR
      INTEGER :: ADD
      INTEGER :: K, I, II, JJ
      integer :: buff_size
      integer :: IAD2, LGTH, IV, IE
      integer :: num_mixed_cells, icell
      integer :: NODEID(4)
      real(kind=WP) :: ALPHA,ALPHAII,ALPHAJJ
      real(kind=WP) :: tol1,tol2
      real(kind=WP) :: NN(3), ALPHAF
      real(kind=WP) :: GRAD(3), NORM_GRAD, VOL_IE
      real(kind=WP) :: PTS(3,4)
      real(kind=WP), pointer, dimension(:) :: ptr_debug


          integer :: mtn,llt,nel,nft,iad,ity,npt,jale,ismstr,jeul,jtur
          integer :: jthe,jlag,jmult,jhbe,jivf,nvaux,jpor,jcvt,jclose,jplasol
          integer :: irep,iint,igtp,israt,isrot,icsen,isorth,isorthg,ifailure,jsms

      INTEGER :: FACE_TO_NODE_LOCAL_ID(4, 2), NODEID1, NODEID2

        real(kind=WP) :: AREA_CELL, A_TARGET
        real(kind=WP) :: D_OUT, D_VAL
        real(kind=WP) :: Dmin, Dmax
        real(kind=WP) ::  D(4) !< projected points on interface normal
        real(kind=WP) :: D_N1, D_N2, EDGE_LEN, T_CLIP
        real(kind=WP) :: Y1, Z1, Y2, Z2
        real(kind=WP) :: VFRAC(21)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      IF(TRIMAT==TRIMAT)ptr_debug => ALE%VOF%cell_data%ALPHA(1:NUMELQ+NUMELS+NQVOIS+NSVOIS,1)

      FACE_TO_NODE_LOCAL_ID(1, 1) = 1 ; FACE_TO_NODE_LOCAL_ID(1, 2) = 2   !!! Face 1
      FACE_TO_NODE_LOCAL_ID(2, 1) = 2 ; FACE_TO_NODE_LOCAL_ID(2, 2) = 3   !!! Face 2
      FACE_TO_NODE_LOCAL_ID(3, 1) = 3 ; FACE_TO_NODE_LOCAL_ID(3, 2) = 4   !!! Face 3
      FACE_TO_NODE_LOCAL_ID(4, 1) = 4 ; FACE_TO_NODE_LOCAL_ID(4, 2) = 1   !!! Face 4

      ! VOLUME FRACTION FOR EACH CELL
      DO NG=ITASK+1,NGROUP,NTHREAD      
!     ALE ON / OFF
         IF (IPARG(76, NG)  ==  1) CYCLE ! --> OFF
         CALL INITBUF(IPARG    ,NG      ,       &
                      MTN     ,LLT     ,NFT     ,IAD     ,ITY     , &
                      NPT     ,JALE    ,ISMSTR  ,JEUL    ,JTUR    , &
                      JTHE    ,JLAG    ,JMULT   ,JHBE    ,JIVF    , &
                      NVAUX   ,JPOR    ,JCVT    ,JCLOSE  ,JPLASOL , &
                      IREP    ,IINT    ,IGTP   ,ISRAT   ,ISROT   ,  &
                      ICSEN   ,ISORTH  ,ISORTHG ,IFAILURE,JSMS    )
         IF(JALE+JEUL == 0) CYCLE
         IF(IPARG(8,NG) == 1) CYCLE
         IF(IPARG(1,NG)  /= 51) CYCLE
         VOLG => ELBUF_TAB(NG)%GBUF%VOL
         UVAR => ELBUF_TAB(NG)%BUFLY(1)%MAT(1,1,1)%VAR
         ! ---  PHASES IN EACH CELL (GLOBAL ARRAY FOR GRADIENT CALCULATION) ---
         DO ITRIMAT = 1, TRIMAT
            ADD    = M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS ! ADD => SIG(1)
            ADD    = ADD + 11   ! ADD + 11 => VOLUME_Phase
            K      = LLT*(ADD-1) ! VAR(I,ADD) = VAR(K+I) 
            VOLP   =>UVAR(K+1:K+LLT)
            !!!   Volume fraction
            DO I=1,LLT
               II     = I+NFT
               ALPHA = VOLP(I)/VOLG(I)
               ALPHA = MAX(ZERO,MIN(ONE,ALPHA))
               ALE%VOF%cell_data%ALPHA(II,ITRIMAT) = ALPHA
               ALE%VOF%cell_data%VOL(II) = VOLG(I)
            ENDDO  
         ENDDO

         ! --- IDENTIFY PRINCIPAL PHASE IN EACH CELL (PHASE WITH LARGEST VOLUME FRACTION) ---
         DO I=1,LLT
           II     = I+NFT
           VFRAC(1:TRIMAT) =  ALE%VOF%cell_data%ALPHA(II,1:TRIMAT)
           ALE%VOF%cell_data%IPHASE(II) = MAXLOC(VFRAC(1:TRIMAT),1)
         ENDDO

      ENDDO  ! NG=ITASK+1,NGROUP,NTHREAD

      CALL MY_BARRIER
      
      !!! MPI Comm
      ! --- EXCHANCE PHASE AT DOMAIN BOUNDARIES
      IF(NSPMD > 1)THEN
        if(itask==0) call startime(timers,timer_spmdcfd)
!$OMP SINGLE
         !!! Volumic fractions comm
         DO ITRIMAT = 1, TRIMAT
            CALL SPMD_E1VOIS(ALE%VOF%cell_data%ALPHA(1,ITRIMAT), NERCVOIS, NESDVOIS,LERCVOIS, LESDVOIS, LENCOM)
         ENDDO
!$OMP END SINGLE
         if(itask==0) call stoptime(timers,timer_spmdcfd)
      ENDIF
      CALL MY_BARRIER 


     ! --- IDENTIFY MIXED CELLS & SET ALPHA AT FACES
     ! A tolerance is introduced, otherwise random PLIC interface => spurious / cell detachment.
     num_mixed_cells = 0
     tol1=1.0e-8_WP
     tol2=ONE - tol1
     ALE%VOF%cell_data%mixed_cell_id(1:NUMELQ+NUMELS) = 0
      DO IE=1,NUMELQ+NUMELS
         IF(ALE%VOF%cell_data%VOL(IE) <= ZERO) CYCLE  ! skip non-ALE elements (VOL not set)
         ITRIMAT = ALE%VOF%cell_data%iphase(IE)
         IF(ITRIMAT < 1 .OR. ITRIMAT > TRIMAT) CYCLE  ! skip uninitialised elements
         ALPHAII = ALE%VOF%cell_data%ALPHA(IE,ITRIMAT)
         IF(ALPHAII > tol1 .AND. ALPHAII < tol2)then
           num_mixed_cells = num_mixed_cells + 1
           ALE%VOF%cell_data%elem_id(num_mixed_cells) = IE
           ALE%VOF%cell_data%mixed_cell_id(IE) = num_mixed_cells
           IAD2 = ALE_CONNECT%ee_connect%iad_connect(IE)
           LGTH = ALE_CONNECT%ee_connect%iad_connect(IE+1) - IAD2
           DO JJ=1,LGTH
              IV = ALE_CONNECT%ee_connect%connected(IAD2 + JJ - 1)
              IF(IV == 0) THEN
                ALPHAJJ = ALPHAII
              ELSEIF(IV > 0) THEN
                ALPHAJJ = ALE%VOF%cell_data%ALPHA(IV,ITRIMAT)
              ELSE
                ALPHAJJ = SEGVAR%PHASE_ALPHA( ITRIMAT,-IV)
              ENDIF
              ALPHA = HALF*(ALPHAII + ALPHAJJ)
              ALE%VOF%cell_data%ALPHA_F(JJ, num_mixed_cells) = ALPHA
           END DO
         endif
     ENDDO
     ALE%VOF%cell_data%num_mixed = num_mixed_cells


     ! --- GRADIENT CALCULATION
      DO ICELL=1,num_mixed_cells
         GRAD(1:3) = ZERO
         IE = ALE%VOF%cell_data%elem_id(ICELL)
         DO JJ=1,4
            NODEID1 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(JJ, 1), IE)
            NODEID2 = IXQ(1 + FACE_TO_NODE_LOCAL_ID(JJ, 2), IE)
            !! compute outward normal vector for current edge. elem connectivity 1,2,3,4 in plane y,z is such as elem normal is in x direction
            !! NODEID1 coordinates are X(1:3, NODEID1) and NODEID2 coordinates are X(1:3, NODEID2)
            !! x is always 0 in 2D
            !! outward normal is NN(2:3)
            NN(1:3) = ZERO
            NN(2) = X(3, NODEID2) - X(3, NODEID1)
            NN(3) = X(2, NODEID1) - X(2, NODEID2)
            ALPHAF = ALE%VOF%cell_data%ALPHA_F(JJ, ICELL)
            GRAD(2) = GRAD(2) + NN(2) * ALPHAF
            GRAD(3) = GRAD(3) + NN(3) * ALPHAF
         END DO
         !! normalize gradient to get unit interface normal direction
         !! GRAD = Sum_f (alpha_f * N_f) where N_f is the outward edge normal (not unit)
         !! The interface normal is n = GRAD / |GRAD| (volume factor cancels out)
         NORM_GRAD = SQRT(GRAD(2)**2 + GRAD(3)**2)
         IF(NORM_GRAD > 1.0e-10_WP) THEN
           ALE%VOF%cell_data%n(1, ICELL) = ZERO
           ALE%VOF%cell_data%n(2:3, ICELL) = GRAD(2:3) / NORM_GRAD
         ELSE
           ! Gradient is negligible (uniform alpha neighborhood) → mark with zero normal.
           ALE%VOF%cell_data%n(1:3, ICELL) = ZERO
         ENDIF
     ENDDO


     !--- PLIC RECONSTRUCTION
     ! FIND D VALUE SUCH AS CLIPPED_AREA_QUAD4(PTS, NN, d) = alpha * AREA_CELL
     ! Using Brent's method on f(d) = CLIPPED_AREA_QUAD4(d) - A_TARGET
     ! f(Dmin) = AREA_CELL - A_TARGET > 0  (full cell at Dmin)
     ! f(Dmax) = 0         - A_TARGET < 0  (empty cell at Dmax)
     ! f is continuous and monotonically decreasing → unique root
     PTS(1, 1:4) = ZERO ! y,z only in 2d
     NN(1) = ZERO ! y,z only in 2d
     DO ICELL=1,num_mixed_cells
        IE = ALE%VOF%cell_data%elem_id(ICELL)
        NN(2:3) = ALE%VOF%cell_data%n(2:3, ICELL)
        ! Skip PLIC reconstruction if normal is zero (degenerate gradient)
        IF(ABS(NN(2)) + ABS(NN(3)) < 1.0e-10_WP) THEN
          ALE%VOF%cell_data%d(ICELL) = ZERO
          CYCLE
        ENDIF
        NODEID(1:4) = IXQ(2:5, IE)
        PTS(2, 1:4) = X(2, NODEID(1:4))
        PTS(3, 1:4) = X(3, NODEID(1:4))
        D(1) = NN(2) * PTS(2, 1) + NN(3) * PTS(3, 1)
        D(2) = NN(2) * PTS(2, 2) + NN(3) * PTS(3, 2)
        D(3) = NN(2) * PTS(2, 3) + NN(3) * PTS(3, 3)
        D(4) = NN(2) * PTS(2, 4) + NN(3) * PTS(3, 4)
        Dmin = MIN(D(1), D(2), D(3), D(4))
        Dmax = MAX(D(1), D(2), D(3), D(4))
        ITRIMAT = ALE%VOF%cell_data%IPHASE(IE)
        ALPHA     = ALE%VOF%cell_data%ALPHA(IE,ITRIMAT)
        AREA_CELL = ALE%VOF%cell_data%VOL(IE)
        A_TARGET  = ALPHA * AREA_CELL
        CALL BRENT_FIND_D(PTS, NN, Dmin, Dmax, A_TARGET, D_OUT)
        ALE%VOF%cell_data%d(ICELL) = D_OUT
     ENDDO


     END SUBROUTINE ALE51_VOF_RECONSTRUCTION2



! ======================================================================================================================
!! \brief Compute the area of the clipped polygon { x in QUAD | n.x >= d } for a 2D QUAD4 cell
!! \details Uses Sutherland-Hodgman clipping of the quad by the half-plane n.x >= d.
!!          The quad is defined in the (y,z) plane (index 2,3 of PTS).
!!          The clipping produces a convex polygon with 0 to 4 vertices.
!!          Area is computed using the shoelace formula.
!! \param[in]  PTS  (3,4) coordinates of the 4 quad nodes (index 1 unused, 2=y, 3=z)
!! \param[in]  NN   (3)   interface normal vector (index 1 unused, 2=ny, 3=nz)
!! \param[in]  D_CUT clipping distance: retain the half-plane n.x >= D_CUT
!! \param[out] AREA  area of the clipped polygon
! ======================================================================================================================
      SUBROUTINE CLIPPED_AREA_QUAD4(PTS, NN, D_CUT, AREA)
        use precision_mod, only : WP
        use constant_mod,  only : ZERO, HALF
        implicit none
        ! ------------------------------------------------------------------
        !                         dummy arguments
        ! ------------------------------------------------------------------
        real(kind=WP), intent(in)  :: PTS(3, 4)
        real(kind=WP), intent(in)  :: NN(3)
        real(kind=WP), intent(in)  :: D_CUT
        real(kind=WP), intent(out) :: AREA
        ! ------------------------------------------------------------------
        !                         local variables
        ! ------------------------------------------------------------------
        integer, parameter :: NMAX = 5  ! max vertices after clipping a quad by 1 line
        real(kind=WP) :: POLY_Y(NMAX), POLY_Z(NMAX)  ! clipped polygon vertices
        integer       :: NPOLY  ! number of vertices in clipped polygon
        real(kind=WP) :: D_I, D_J, T
        real(kind=WP) :: YI, ZI, YJ, ZJ
        integer       :: II, JJ
        real(kind=WP) :: SUM_CROSS
        ! ------------------------------------------------------------------
        !                         body
        ! ------------------------------------------------------------------

        ! --- Sutherland-Hodgman clipping: quad (4 edges) vs half-plane n.x >= D_CUT ---
        ! Input polygon = quad with 4 vertices
        ! Output polygon = clipped polygon with 0..5 vertices
        NPOLY = 0
        DO II = 1, 4
          JJ = MOD(II, 4) + 1  ! next vertex (cyclic: 1->2->3->4->1)
          YI = PTS(2, II) ; ZI = PTS(3, II)
          YJ = PTS(2, JJ) ; ZJ = PTS(3, JJ)
          D_I = NN(2) * YI + NN(3) * ZI - D_CUT  ! signed distance of vertex I to clipping line
          D_J = NN(2) * YJ + NN(3) * ZJ - D_CUT  ! signed distance of vertex J to clipping line

          IF(D_I >= ZERO) THEN
            ! vertex I is INSIDE (n.x >= d)
            NPOLY = NPOLY + 1
            POLY_Y(NPOLY) = YI
            POLY_Z(NPOLY) = ZI
            IF(D_J < ZERO) THEN
              ! edge I->J crosses the line (I inside, J outside) → add intersection
              T = D_I / (D_I - D_J)
              NPOLY = NPOLY + 1
              POLY_Y(NPOLY) = YI + T * (YJ - YI)
              POLY_Z(NPOLY) = ZI + T * (ZJ - ZI)
            ENDIF
          ELSE
            ! vertex I is OUTSIDE (n.x < d)
            IF(D_J >= ZERO) THEN
              ! edge I->J crosses the line (I outside, J inside) → add intersection
              T = D_I / (D_I - D_J)
              NPOLY = NPOLY + 1
              POLY_Y(NPOLY) = YI + T * (YJ - YI)
              POLY_Z(NPOLY) = ZI + T * (ZJ - ZI)
            ENDIF
          ENDIF
        ENDDO

        ! --- Compute area of clipped polygon using shoelace formula ---
        IF(NPOLY < 3) THEN
          AREA = ZERO
          RETURN
        ENDIF

        SUM_CROSS = ZERO
        DO II = 1, NPOLY
          JJ = MOD(II, NPOLY) + 1
          SUM_CROSS = SUM_CROSS + POLY_Y(II) * POLY_Z(JJ) - POLY_Y(JJ) * POLY_Z(II)
        ENDDO
        AREA = HALF * ABS(SUM_CROSS)

      END SUBROUTINE CLIPPED_AREA_QUAD4


! ======================================================================================================================
!! \brief Find d in [Dmin, Dmax] such that CLIPPED_AREA_QUAD4(PTS, NN, d) = A_TARGET using Brent's method
!! \details f(d) = CLIPPED_AREA_QUAD4(d) - A_TARGET is continuous, monotonically decreasing.
!!          f(Dmin) > 0 (full cell), f(Dmax) < 0 (empty cell) → guaranteed root in [Dmin, Dmax].
!!          Convergence in ~10 iterations to machine precision.
!! \param[in]  PTS      (3,4) coordinates of the 4 quad nodes
!! \param[in]  NN       (3)   interface normal vector
!! \param[in]  DMIN     lower bound (smallest node projection)
!! \param[in]  DMAX     upper bound (largest node projection)
!! \param[in]  A_TARGET target clipped area = alpha * AREA_CELL
!! \param[out] D_OUT    d value such that CLIPPED_AREA_QUAD4(d) ≈ A_TARGET
! ======================================================================================================================
      SUBROUTINE BRENT_FIND_D(PTS, NN, DMIN, DMAX, A_TARGET, D_OUT)
        use precision_mod, only : WP
        use constant_mod,  only : ZERO, HALF, EM20
        implicit none
        ! ------------------------------------------------------------------
        !                         dummy arguments
        ! ------------------------------------------------------------------
        real(kind=WP), intent(in)  :: PTS(3, 4)
        real(kind=WP), intent(in)  :: NN(3)
        real(kind=WP), intent(in)  :: DMIN, DMAX
        real(kind=WP), intent(in)  :: A_TARGET
        real(kind=WP), intent(out) :: D_OUT
        ! ------------------------------------------------------------------
        !                         local variables
        ! ------------------------------------------------------------------
        integer, parameter :: MAX_ITER = 50
        real(kind=WP), parameter :: TOL_REL = 1.0e-12_WP
        real(kind=WP) :: A, B, C, DD, E, S
        real(kind=WP) :: FA, FB, FC, FS
        real(kind=WP) :: TOL, M
        real(kind=WP) :: AREA_TMP
        integer       :: ITER
        logical       :: USE_BISECT
        ! ------------------------------------------------------------------
        !                         body
        ! ------------------------------------------------------------------

        ! Evaluate f at bounds
        A = DMIN
        CALL CLIPPED_AREA_QUAD4(PTS, NN, A, AREA_TMP)
        FA = AREA_TMP - A_TARGET  ! should be > 0 (full cell - target)

        B = DMAX
        CALL CLIPPED_AREA_QUAD4(PTS, NN, B, AREA_TMP)
        FB = AREA_TMP - A_TARGET  ! should be < 0 (empty cell - target)

        ! Safety: if no sign change, return best bound
        IF(FA * FB > ZERO) THEN
          IF(ABS(FA) < ABS(FB)) THEN
            D_OUT = A
          ELSE
            D_OUT = B
          ENDIF
          RETURN
        ENDIF

        ! Ensure |f(a)| >= |f(b)| (b is the best current guess)
        IF(ABS(FA) < ABS(FB)) THEN
          ! swap a and b
          C = A ; A = B ; B = C
          FC = FA ; FA = FB ; FB = FC
        ENDIF

        C = A ; FC = FA
        DD = B - A
        E  = DD

        DO ITER = 1, MAX_ITER
          ! Convergence check
          TOL = TOL_REL * (ABS(B) + ABS(A) + real(1,WP))
          M   = HALF * (C - B)
          IF(ABS(M) <= TOL .OR. ABS(FB) <= EM20) THEN
            D_OUT = B
            RETURN
          ENDIF

          ! Decide: inverse quadratic interpolation or bisection
          USE_BISECT = .TRUE.
          IF(ABS(E) >= TOL .AND. ABS(FA) > ABS(FB)) THEN
            ! Try inverse quadratic interpolation
            S = FB / FA
            IF(ABS(A - C) <= EM20) THEN
              ! Secant method (a == c)
              DD = real(2,WP) * M * S / (real(1,WP) - S)
            ELSE
              ! Inverse quadratic interpolation
              DD = S * (real(2,WP) * M * FC * (FC - FB) - (B - A) * (FB - FC))
              DD = DD / ((FA - FB) * (FA - FC) * (FB - FC) + EM20)
            ENDIF
            ! Accept interpolation step if it stays within bounds
            IF(ABS(DD) < ABS(HALF * E) .AND. &
               DD > real(-1,WP) * ABS(M) + TOL .AND. &
               DD < ABS(M) - TOL) THEN
              USE_BISECT = .FALSE.
              E = DD
            ENDIF
          ENDIF

          IF(USE_BISECT) THEN
            DD = M
            E  = M
          ENDIF

          ! Update a (previous best)
          A  = B
          FA = FB

          ! Update b (new best)
          IF(ABS(DD) > TOL) THEN
            B = B + DD
          ELSE
            ! minimum step
            IF(M > ZERO) THEN
              B = B + TOL
            ELSE
              B = B - TOL
            ENDIF
          ENDIF

          CALL CLIPPED_AREA_QUAD4(PTS, NN, B, AREA_TMP)
          FB = AREA_TMP - A_TARGET

          ! Maintain bracket: c and b must have opposite signs
          IF(FB * FC > ZERO) THEN
            C  = A  ; FC = FA
            DD = B - A ; E = DD
          ENDIF

          ! Ensure |f(b)| <= |f(c)|
          IF(ABS(FC) < ABS(FB)) THEN
            A = B ; B = C ; C = A
            FA = FB ; FB = FC ; FC = FA
          ENDIF
        ENDDO

        D_OUT = B

      END SUBROUTINE BRENT_FIND_D



     END MODULE ALE51_VOF_RECONSTRUCTION2_MOD

