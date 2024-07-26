Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
Chd|====================================================================
Chd|  EBCS11                        source/boundary_conditions/ebcs/ebcs10.F
Chd|-- called by -----------
Chd|        EBCS_MAIN                     source/boundary_conditions/ebcs/ebcs_main.F
Chd|-- calls ---------------
Chd|        ALE_MOD                       ../common_source/modules/ale/ale_mod.F
Chd|        EBCS_MOD                      ../common_source/modules/boundary_conditions/ebcs_mod.F
Chd|        ELBUFDEF_MOD                  ../common_source/modules/mat_elem/elbufdef_mod.F
Chd|        MULTI_FVM_MOD                 ../common_source/modules/ale/multi_fvm_mod.F
Chd|        SEGVAR_MOD                    share/modules/segvar_mod.F    
Chd|        TH_SURF_MOD                   ../common_source/modules/interfaces/th_surf_mod.F
Chd|====================================================================
      SUBROUTINE EBCS11(NSEG,ISEG,SEGVAR,
     .                  V,W,X,
     .                  LISTE,NOD,IRECT,IELEM,IFACE,
     .                  LA,MS,STIFN,EBCS,IPARG,ELBUF_TAB,IXQ,IXS,IXTG,
     .                  FSAVSURF)
      USE EBCS_MOD
      USE ELBUFDEF_MOD
      USE MULTI_FVM_MOD
      USE SEGVAR_MOD
      USE ALE_MOD , only : ALE
      USE TH_SURF_MOD , only : TH_SURF_NUM_CHANNEL
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "param_c.inc"
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "com08_c.inc"
#include      "mmale51_c.inc"
#include      "parit_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      my_real, INTENT(INOUT) :: FSAVSURF(TH_SURF_NUM_CHANNEL,NSURF)
      INTEGER,INTENT(IN) :: NSEG,NOD,ISEG(NSEG),LISTE(NOD),IRECT(4,NSEG),IELEM(NSEG),IFACE(NSEG)
      INTEGER,INTENT(IN) :: IXQ(NIXQ,NUMELQ),IXS(NIXS,NUMELS),IXTG(NIXTG,NUMELTG)
      my_real,INTENT(INOUT) :: MS(NUMNOD)
      my_real V(3,NUMNOD),W(3,NUMNOD),X(3,NUMNOD),LA(3,NOD),STIFN(NUMNOD)
      TYPE(t_ebcs_propergol), INTENT(INOUT) :: EBCS
      INTEGER :: IPARG(NPARG,NGROUP)
      TYPE(ELBUF_STRUCT_), TARGET, DIMENSION(NGROUP) :: ELBUF_TAB
      TYPE(t_segvar),INTENT(INOUT) :: SEGVAR
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      TYPE(G_BUFEL_), POINTER :: GBUF
      TYPE(L_BUFEL_)  ,POINTER :: LBUF     
      INTEGER II,IS,KK,K1,LL,KSEG,NN(4),NUM,KTY,KLT,MFT,NGRP,ILOC,NPT,IVOI,IDX(6),IX(4)
      INTEGER ICF_2d(2,4), ICF_3d(4,6), JJ, ISUBMAT, IPOS, MTN
      my_real ORIENT,RHO,C,FC,ROC,ALP,VOL,MASS,
     .        X13,Y13,Z13,X24,Y24,Z24,NX,NY,Z,S,
     .        ROOU,ENOU,VMX,VMY,VMZ,FLUXI,FLUXO,P,DVX,DVY,DVZ,
     .        V0(3,NOD),
     .        XSN,YSN,ZSN,
     .        XN, YN, ZN, Vold, Vnew, Pold, Padj,
     .        DP0,MACH,PP,SSP,RHOC2,VEL_FRONT,FAC1,FAC2,
     .        PARAM_A, PARAM_N, PARAM_Q, SURF, TIME, EINT,PHASE_ALPHA(21),PHASE_RHO(21), PHASE_EINT(21)
      my_real :: RHO0S     !< initial propergol density
      my_real :: DMASS_G   !< mass increment (burnt propergol)
      my_real :: DVOL_S    !< burnt volume
      my_real :: DEINT_G   !< released energy from burnt propergol

      LOGICAL bFOUND
      TYPE(BUF_MAT_)  ,POINTER :: MBUF  
      INTEGER :: ADRESS  !< adress for parit/on array
      INTEGER :: SURF_ID !< EBCS surface identifier
      INTEGER :: IMAT    !< submaterial identifier for burnt gas
C-----------------------------------------------
      DATA ICF_2d  /1,2,2,3,3,4,4,1/
      DATA ICF_3d  /1,4,3,2,3,4,8,7,5,6,7,8,1,2,6,5,2,3,7,6,1,5,8,4/
C-----------------------------------------------
C   S o u r  c e   L i n e s
C-----------------------------------------------
      PARAM_A = EBCS%A
      PARAM_N = EBCS%N
      PARAM_Q = EBCS%Q
      TIME = TT
      SURF_ID = EBCS%SURF_ID

      !vitesse relative et modele de combustion : nettoyage possible ?

      DO II=1,NOD
        NUM=LISTE(II)
        V0(1,II)=V(1,NUM)
        V0(2,II)=V(2,NUM)
        V0(3,II)=V(3,NUM)
      ENDDO

      IF(IALE == 1)THEN
        DO II=1,NOD
          NUM=LISTE(II)
          V0(1,II)=V0(1,II)-W(1,NUM)
          V0(2,II)=V0(1,II)-W(2,NUM)
          V0(3,II)=V0(1,II)-W(3,NUM)
        ENDDO
      ENDIF

      DO II=1,NOD
        NUM=LISTE(II)
        LA(1,II)=ZERO
        LA(2,II)=ZERO
        LA(3,II)=ZERO
      ENDDO

        DO IS=1,NSEG
        
          KSEG=ABS(ISEG(IS))
          ORIENT=FLOAT(ISEG(IS)/KSEG)
                   
          !---OUTWARD NORMAL
          IF(N2D == 0)THEN
            JJ = IFACE(IS)
            IX(1)=IXS(ICF_3d(1,JJ)+1,IELEM(IS))          
            IX(2)=IXS(ICF_3d(2,JJ)+1,IELEM(IS))          
            IX(3)=IXS(ICF_3d(3,JJ)+1,IELEM(IS))          
            IX(4)=IXS(ICF_3d(4,JJ)+1,IELEM(IS)) 
            X13=X(1,IX(3))-X(1,IX(1))           
            Y13=X(2,IX(3))-X(2,IX(1))           
            Z13=X(3,IX(3))-X(3,IX(1))           
            X24=X(1,IX(4))-X(1,IX(2))           
            Y24=X(2,IX(4))-X(2,IX(2))           
            Z24=X(3,IX(4))-X(3,IX(2))           
            XN=Y13*Z24-Z13*Y24             
            YN=Z13*X24-X13*Z24             
            ZN=X13*Y24-Y13*X24             
            FAC2=ONE/SQRT(XN**2+YN**2+ZN**2)
            XN = XN*FAC2
            YN = YN*FAC2
            ZN = ZN*FAC2
            SURF = HALF/FAC2
            IF(IX(4) == IX(3))THEN ; NPT=3;FAC1=THIRD; ELSE; NPT=4;FAC1=FOURTH; ENDIF
          ELSE
            FAC1=HALF
            NPT=2
            JJ = IFACE(IS)
            IF(NUMELTG>0)THEN
              IX(1)  = IXTG(ICF_2d(1,JJ)+1,IELEM(IS))
              IX(2)  = IXTG(ICF_2d(2,JJ)+1,IELEM(IS))
            ELSE
              IX(1)  = IXQ(ICF_2d(1,JJ)+1,IELEM(IS))
              IX(2)  = IXQ(ICF_2d(2,JJ)+1,IELEM(IS))
            ENDIF 
            XN = ZERO
            YN = -(-X(3,IX(2))+X(3,IX(1)))
            ZN =  (-X(2,IX(2))+X(2,IX(1)))                     
            FAC2 = ONE/SQRT(YN*YN+ZN*ZN)
            YN=YN*FAC2
            ZN=ZN*FAC2
            SURF = ONE/FAC2
          ENDIF

          !-- ADJACENT STATE
          bFOUND = .FALSE.
          IVOI = IELEM(IS)                                  
          DO NGRP=1,NGROUP                                     
            KTY = IPARG(5,NGRP)                                
            KLT = IPARG(2,NGRP)                                
            MFT = IPARG(3,NGRP)
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

          !PRESSURE                                         
          ILOC = IVOI-MFT-1
          DO KK=1,6                                          
            IDX(KK) = KLT*(KK-1)                               
          ENDDO                                             
          Padj = -THIRD*(GBUF%SIG(IDX(1)+ILOC+1) + GBUF%SIG(IDX(2)+ILOC+1) + GBUF%SIG(IDX(3)+ILOC+1))

          !DENSITY                                          
          RHO = GBUF%RHO(ILOC+1)

          !VOLUME
          VOL = GBUF%VOL(ILOC+1)
          MASS = RHO*VOL

          !ENERGY
          EINT = GBUF%EINT(ILOC+1)

          !SOUND SPEED
          SSP = LBUF%SSP(ILOC+1)
          
          !VOL FRAC          
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
          ENDIF



          !-- FORMULATION

          !FRONT VELOCITY
          PP = Padj
          IF(PP <= ZERO)THEN
            VEL_FRONT = ZERO
            DVOL_S = ZERO
          ELSE
            VEL_FRONT = PARAM_A*EXP(PARAM_N*LOG(PP))
            DVOL_S = SURF*VEL_FRONT*DT1
          ENDIF

          !INCREMENT (PROPERGOL -> GAS)
          DMASS_G = RHO0S * DVOL_S
          DEINT_G = PARAM_Q * DMASS_G

          ! element state
          RHO = (RHO*VOL + DMASS_G) / VOL !(VOL + DVOL_S)
          EINT = (EINT*VOL + DEINT_G) / VOL
          ! element update
          GBUF%RHO(ILOC+1) = RHO
          GBUF%EINT(ILOC+1) = EINT

          ! burnt gas is supposed to be submat 1
          IMAT = 1
          IF(SEGVAR%has_phase_rho) THEN
            VOL = PHASE_ALPHA(IMAT)*VOL
            PHASE_RHO(IMAT) = (PHASE_RHO(IMAT)*VOL + DMASS_G) / VOL
            PHASE_EINT(IMAT) = (PHASE_EINT(IMAT)*VOL + DEINT_G) / VOL
          ENDIF

          ! ghost cell update ( upwind/ACONVE() )
          SEGVAR%RHO(KSEG)=RHO
          SEGVAR%EINT(KSEG)=EINT
          IF(SEGVAR%has_phase_alpha)SEGVAR%PHASE_ALPHA(1:4,KSEG)=PHASE_ALPHA(1:4)
          IF(SEGVAR%has_phase_rho)SEGVAR%PHASE_RHO(1:4,KSEG)=PHASE_RHO(1:4)
          IF(SEGVAR%has_phase_eint)SEGVAR%PHASE_EINT(1:4,KSEG)=PHASE_EINT(1:4)

          !-- NODES
          NN(1)=IRECT(1,IS)
          NN(2)=IRECT(2,IS)
          NN(3)=IRECT(3,IS)
          NN(4)=IRECT(4,IS)

          !expand pressure loading to segment nodes
          DO KK=1,NPT

            LA(1,NN(KK)) = LA(1,NN(KK)) -  PP*SURF*XN*FAC1
            LA(2,NN(KK)) = LA(2,NN(KK)) -  PP*SURF*YN*FAC1
            LA(3,NN(KK)) = LA(3,NN(KK)) -  PP*SURF*ZN*FAC1

            IF (MS(NN(KK)) > ZERO)THEN
              STIFN(NN(KK))=STIFN(NN(KK))+(TWO*(SURF*RHO*SSP)**2)/MS(NN(KK))
            ENDIF

          ENDDO

           !/TH/SURF (massflow, velocity, pressure)
           FSAVSURF(2,SURF_ID) = FSAVSURF(2,SURF_ID) + RHO*SURF*VEL_FRONT     !rho.S.u = dm/dt
           FSAVSURF(3,SURF_ID) = FSAVSURF(3,SURF_ID) + SURF*VEL_FRONT         !S.u
           FSAVSURF(4,SURF_ID) = FSAVSURF(4,SURF_ID) + SURF*Padj              !S.P
           FSAVSURF(6,SURF_ID) = FSAVSURF(6,SURF_ID) + RHO*SURF*VEL_FRONT*DT1 ! m<-m+dm (cumulative)

        ENDDO      
        
        !no force resulting from ghost cell : Pelem = Padj (burnt gas), Pghost=0 (inert solid)



      RETURN
      END
