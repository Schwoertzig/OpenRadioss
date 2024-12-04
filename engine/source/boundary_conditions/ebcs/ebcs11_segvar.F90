!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
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
!                                                   PROCEDURES
! ======================================================================================================================
      SUBROUTINE EBCS11_SEGVAR(NSEG,ISEG,SEGVAR,NOD,IELEM,EBCS,IPARG,ELBUF_TAB,DT1,NPARG,NGROUP,N2D)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE EBCS_MOD
      USE ELBUFDEF_MOD
      USE MULTI_FVM_MOD
      USE SEGVAR_MOD
      USE TH_SURF_MOD , only : TH_SURF_NUM_CHANNEL
      USE CONSTANT_MOD , only : ZERO, THIRD
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit None
! ----------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER, INTENT(IN) :: N2D
      INTEGER,INTENT(IN) :: NPARG, NGROUP !< array sizes
      my_real, INTENT(IN) :: DT1 !< time step
      INTEGER,INTENT(IN) :: NSEG,NOD,ISEG(NSEG),IELEM(NSEG)
      TYPE(t_ebcs_propergol), INTENT(INOUT) :: EBCS
      INTEGER :: IPARG(NPARG,NGROUP)
      TYPE(ELBUF_STRUCT_), TARGET, DIMENSION(NGROUP) :: ELBUF_TAB
      TYPE(t_segvar),INTENT(INOUT) :: SEGVAR

      INTEGER N0PHAS
      PARAMETER (N0PHAS = 04)

      INTEGER NVPHAS
      PARAMETER (NVPHAS = 23)

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      TYPE(G_BUFEL_), POINTER :: GBUF
      TYPE(L_BUFEL_),POINTER :: LBUF
      INTEGER IS,KK,KSEG,NN(4),NNG(4),NUM,KTY,KLT,MFT,NGRP,ILOC,NPT,IVOI,IDX(6)
      INTEGER ICF_2d(2,4), ICF_3d(4,6), ISUBMAT, IPOS, MTN
      INTEGER ISOLNOD
      my_real :: RHO,VOL,MASS,MASS_FACE,V0(3,NOD),Vold, Vnew, Pold, Padj,PP,SSP,RHOC2,VEL_FRONT
      my_real :: PARAM_A, PARAM_N, PARAM_Q, SURF,PHASE_ALPHA(21),PHASE_RHO(21), PHASE_EINT(21)
      my_real :: PARAM_RHO0S  !< initial propergol density
      my_real :: DMASS_G     !< mass increment (burnt propergol)
      my_real :: DVOL_S      !< burnt volume
      my_real :: DEINT_G     !< released energy from burnt propergol
      my_real :: TMP(3)

      LOGICAL bFOUND
      TYPE(BUF_MAT_)  ,POINTER :: MBUF
      INTEGER :: SURF_ID !< EBCS surface identifier
      INTEGER :: NBSUBMAT
      !-----------------------------------------------
      DATA ICF_2d  /1,2,2,3,3,4,4,1/
      DATA ICF_3d  /1,4,3,2,3,4,8,7,5,6,7,8,1,2,6,5,2,3,7,6,1,5,8,4/
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      PARAM_A = EBCS%A
      PARAM_N = EBCS%N
      PARAM_Q = EBCS%Q
      PARAM_RHO0S = EBCS%rho0s
      SURF_ID = EBCS%SURF_ID


        DO IS=1,NSEG
        
          KSEG=ABS(ISEG(IS))

          !-- ADJACENT STATE
          bFOUND = .FALSE.
          IVOI = IELEM(IS)                                  
          DO NGRP=1,NGROUP                                     
            KTY = IPARG(5,NGRP)                                
            KLT = IPARG(2,NGRP)                                
            MFT = IPARG(3,NGRP)
            ISOLNOD = IPARG(28,NGRP)
            IF(N2D==0)THEN
              IF(KTY /= 1)CYCLE
            ELSE
              IF(KTY /= 2 .AND. KTY /= 7)CYCLE
            ENDIF                                
             IF (IVOI <= KLT+MFT)THEN        
               bFOUND = .TRUE.                              
               EXIT                                         
             ENDIF                                          
          ENDDO
          IF(.NOT.bFOUND)CYCLE !next I                 
          GBUF => ELBUF_TAB(NGRP)%GBUF
          LBUF => ELBUF_TAB(NGRP)%BUFLY(1)%LBUF(1,1,1)
          MTN = IPARG(1,NGRP)
          !ADJACENT PRESSURE
          ILOC = IVOI-MFT-1
          DO KK=1,6                                          
            IDX(KK) = KLT*(KK-1)                               
          ENDDO                                             
          Padj = -THIRD*(GBUF%SIG(IDX(1)+ILOC+1) + GBUF%SIG(IDX(2)+ILOC+1) + GBUF%SIG(IDX(3)+ILOC+1))
          !DENSITY
          RHO = GBUF%RHO(ILOC+1)
          VOL = GBUF%VOL(ILOC+1)
          MASS = VOL * RHO
          !ADJACENT SOUND SPEED
          SSP = LBUF%SSP(ILOC+1)
          !VOLUME FRACIONS AND SUBMAT STATE
          IF(MTN == 51)THEN
            MBUF => ELBUF_TAB(NGRP)%BUFLY(1)%MAT(1,1,1)
            DO ISUBMAT=1,4
              IPOS = 1                                                                                                       
              KK = (N0PHAS + (ISUBMAT-1)*NVPHAS +IPOS-1) * KLT  +  ILOC+1                                                                                    
              PHASE_ALPHA(ISUBMAT) = MBUF%VAR(KK)                                     
              IPOS = 9                                                                                                       
              KK = (N0PHAS + (ISUBMAT-1)*NVPHAS +IPOS-1) * KLT  +  ILOC+1                                                                                    
              PHASE_RHO(ISUBMAT) = MBUF%VAR(KK)
              IPOS = 8                                                                                                       
              KK = (N0PHAS + (ISUBMAT-1)*NVPHAS +IPOS-1) * KLT  +  ILOC+1                                                                                    
              PHASE_EINT(ISUBMAT) = MBUF%VAR(KK)
            ENDDO!next ITRIMAT
            VOL = VOL * PHASE_ALPHA(1) ! burnt gas is supposed to be submat 1
            MASS = VOL * PHASE_RHO(1)
          ENDIF

          !-- FORMULATION
          !BURNT PROPERGOL VOLUME
          PP = Padj
           !adjacent elem
          IF(PP <= ZERO)THEN
            VEL_FRONT = ZERO
            DVOL_S = ZERO
            DMASS_G = ZERO
          ELSE
            VEL_FRONT = PARAM_A*EXP(PARAM_N*LOG(PP))
            DVOL_S = SURF*VEL_FRONT*DT1
            DMASS_G = DVOL_S * PARAM_RHO0S
          ENDIF


          ! ghost cell update ( upwind/ACONVE() )
          SEGVAR%RHO(KSEG)=PARAM_RHO0S
          SEGVAR%EINT(KSEG)=PARAM_Q*PARAM_RHO0S


          ! burnt gas is supposed to be submat 1
          IF(SEGVAR%nbmat > 1)THEN
            NBSUBMAT = SEGVAR%NBMAT

            PHASE_RHO(1) = PARAM_RHO0S
            SEGVAR%PHASE_RHO(2:NBSUBMAT,KSEG)=PHASE_RHO(2:NBSUBMAT)

            PHASE_EINT(1) = PARAM_Q*PARAM_RHO0S
            SEGVAR%PHASE_EINT(2:NBSUBMAT,KSEG)=PHASE_EINT(2:NBSUBMAT)

            SEGVAR%PHASE_ALPHA(1:NBSUBMAT,KSEG)=PHASE_ALPHA(1:NBSUBMAT)
          ENDIF

          ! -------------
        ENDDO




      RETURN
      END
