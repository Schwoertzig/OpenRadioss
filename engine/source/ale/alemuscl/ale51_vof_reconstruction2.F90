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
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
      SUBROUTINE ALE51_VOF_RECONSTRUCTION2(IPARG   ,ELBUF_TAB      ,ALE_CONNECT, &
                                           NERCVOIS ,NESDVOIS,LERCVOIS,LESDVOIS   ,LENCOM, ITASK, &
                                           SEGVAR  ,timers , &
                                           NPARG,NGROUP,NUMELQ,NUMELS,TRIMAT,NQVOIS,NSVOIS,NSPMD)
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
      use constant_mod , only : ZERO, ONE, EP20, FOURTH, HALF
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
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: NG
      INTEGER :: ITRIMAT
      real(kind=WP), DIMENSION(:), POINTER :: VOLG, VOLP, UVAR
      INTEGER :: ADD
      INTEGER :: K, I, II, JJ
      INTEGER :: ELEM_ID
      integer :: buff_size
      integer :: IAD2, LGTH, IV, IE
      real(kind=WP) :: ALPHA,ALPHAII,ALPHAJJ
      real(kind=WP), pointer, dimension(:) :: ptr_debug


          integer :: mtn,llt,nel,nft,iad,ity,npt,jale,ismstr,jeul,jtur
          integer :: jthe,jlag,jmult,jhbe,jivf,nvaux,jpor,jcvt,jclose,jplasol
          integer :: irep,iint,igtp,israt,isrot,icsen,isorth,isorthg,ifailure,jsms
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      IF(TRIMAT==TRIMAT)ptr_debug => ALE%VOF%cell_data%ALPHA(1:NUMELQ+NUMELS+NQVOIS+NSVOIS)


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
         DO ITRIMAT = 1, 1
            ADD    = M51_N0PHAS + (ITRIMAT-1)*M51_NVPHAS ! ADD => SIG(1)
            ADD    = ADD + 11   ! ADD + 11 => VOLUME_Phase
            K      = LLT*(ADD-1) ! VAR(I,ADD) = VAR(K+I) 
            VOLP   =>UVAR(K+1:K+LLT)
            !!!   Volume fraction
            DO I=1,LLT
               II     = I+NFT
               ALPHA = VOLP(I)/VOLG(I)
               print *, "---" , ALPHA
               ALPHA = MAX(ZERO,MIN(ONE,ALPHA))
               ALE%VOF%cell_data%ALPHA(II) = ALPHA
            ENDDO  
         ENDDO  
      ENDDO  ! NG=ITASK+1,NGROUP,NTHREAD

      CALL MY_BARRIER
      
      !!! MPI Comm
      IF(NSPMD > 1)THEN
        if(itask==0) call startime(timers,timer_spmdcfd)
!$OMP SINGLE
         !!! Volumic fractions comm
         DO ITRIMAT = 1, 1
            CALL SPMD_E1VOIS(ALE%VOF%cell_data%ALPHA(1), NERCVOIS, NESDVOIS,LERCVOIS, LESDVOIS, LENCOM)
         ENDDO
!$OMP END SINGLE
         if(itask==0) call stoptime(timers,timer_spmdcfd)
      ENDIF
      CALL MY_BARRIER 

     ITRIMAT = 1
      DO IE=1,NUMELQ+NUMELS
         ALPHAII = ALE%VOF%cell_data%ALPHA(IE)
         IAD2 = ALE_CONNECT%ee_connect%iad_connect(IE)
         LGTH = ALE_CONNECT%ee_connect%iad_connect(IE+1) - IAD2
         DO JJ=1,LGTH
            IV = ALE_CONNECT%ee_connect%connected(IAD2 + JJ - 1)
            IF(IV == 0) THEN
              ALPHAJJ = ALPHAII
            ELSEIF(IV > 0) THEN
              ALPHAJJ = ALE%VOF%cell_data%ALPHA(IV)
            ELSE
              ALPHAJJ = SEGVAR%PHASE_ALPHA( ITRIMAT,-IV)
            ENDIF
            ALPHA = HALF*(ALPHAII + ALPHAJJ)
            ALE%VOF%cell_data%ALPHA_F(JJ, IE) = ALPHA
         END DO
     ENDDO




     END SUBROUTINE ALE51_VOF_RECONSTRUCTION2
     END MODULE ALE51_VOF_RECONSTRUCTION2_MOD
      
