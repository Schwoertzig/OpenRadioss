Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
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
      !||====================================================================
      !||    sigeps109             ../engine/source/materials/mat/mat109/sigeps109.F
      !||--- called by ------------------------------------------------------
      !||    mulaw                 ../engine/source/materials/mat_share/mulaw.F90
      !||--- calls      -----------------------------------------------------
      !||    table2d_vinterp_log   ../engine/source/tools/curve/table2d_vinterp_log.F
      !||    table_vinterp         ../engine/source/tools/curve/table_tools.F
      !||--- uses       -----------------------------------------------------
      !||    interface_table_mod   ../engine/share/modules/table_mod.F
      !||    table_mod             ../engine/share/modules/table_mod.F
      !||====================================================================
      SUBROUTINE SIGEPS109(
     1     NEL     ,NGL     ,NUPARAM ,NUVAR   ,NVARTMP ,NUMTABL ,
     2     UPARAM  ,UVAR    ,VARTMP  ,ITABLE  ,TABLE   ,JTHE    ,
     3     TIME    ,TIMESTEP,OFF     ,RHO     ,PLA     ,DPLA    ,
     4     SOUNDSP ,SIGY    ,ET      ,TEMP    ,EPSP    ,DPDM    ,
     5     DEPSXX  ,DEPSYY  ,DEPSZZ  ,DEPSXY  ,DEPSYZ  ,DEPSZX  ,
     6     SIGOXX  ,SIGOYY  ,SIGOZZ  ,SIGOXY  ,SIGOYZ  ,SIGOZX  ,
     7     SIGNXX  ,SIGNYY  ,SIGNZZ  ,SIGNXY  ,SIGNYZ  ,SIGNZX  ,
     8     INLOC   ,DPLANL  ,IEOS    ,JLAG    ,FHEAT   ,VOLUME  )
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE TABLE_MOD
      USE INTERFACE_TABLE_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
C-----------------------------------------------
C   D u m m y A r g u m e n t s
C-----------------------------------------------
      INTEGER ,INTENT(IN) :: IEOS
      INTEGER ,INTENT(IN) :: JTHE
      INTEGER ,INTENT(IN) :: JLAG
      INTEGER NEL,NUPARAM,NUVAR,NVARTMP,NUMTABL,INLOC
      INTEGER ,DIMENSION(NUMTABL),INTENT(IN)  :: ITABLE
      INTEGER ,DIMENSION(NEL)    ,INTENT(IN)  :: NGL
c
      my_real  :: TIME,TIMESTEP
      my_real,DIMENSION(NUPARAM) ,INTENT(IN)  :: UPARAM
      my_real,DIMENSION(NEL)     ,INTENT(IN)  :: VOLUME
      my_real,DIMENSION(NEL)     ,INTENT(IN)  :: RHO,DPDM,OFF,
     .   DEPSXX,DEPSYY,DEPSZZ,DEPSXY,DEPSYZ,DEPSZX,
     .   SIGOXX,SIGOYY,SIGOZZ,SIGOXY,SIGOYZ,SIGOZX,DPLANL
      my_real ,DIMENSION(NEL)    ,INTENT(OUT) ::  SOUNDSP,SIGY,ET,
     .   SIGNXX,SIGNYY,SIGNZZ,SIGNXY,SIGNYZ,SIGNZX
      my_real ,DIMENSION(NEL) ,INTENT(INOUT) :: PLA,DPLA,EPSP
      my_real ,DIMENSION(NEL,NUVAR)   ,INTENT(INOUT) :: UVAR
      INTEGER ,DIMENSION(NEL,NVARTMP) ,INTENT(INOUT) :: VARTMP
      my_real ,DIMENSION(NEL) ,INTENT(INOUT) :: TEMP
      my_real ,DIMENSION(NEL) ,INTENT(INOUT) :: FHEAT
c
      TYPE(TTABLE), DIMENSION(NTABLE) ::  TABLE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I,II,NINDX,ITER,NITER,ISMOOTH,
     .   FUNC_YLD,FUNC_TEMP,FUNC_ETA,NDIM_YLD,NDIM_TEMP,NDIM_ETA
      INTEGER ,DIMENSION(NEL) :: INDEX
c
      my_real  :: YOUNG,BULK,LAME,G,G2,G3,NU,TREF,TINI,ETA,LDAV,CP,
     .   XRATE,XSCALE,YSCALE,DTINV,J2,Q2,DPHI_DLAM,R,DLAMIN,ALPHA,ALPHI
      my_real  :: check_svm
c
      my_real, DIMENSION(NEL) ::SVM,SVMT,YLD,YLD_TREF,YLD_TEMP,
     .   SXX,SYY,SZZ,SXY,SYZ,SZX,SIGM,STXX,STYY,STZZ,STXY,STYZ,STZX,
     .   FACT_ETA,DYDX,HARDP,HARDR,YLD_I,HARDP_I,HARDR_I,DXDYV,DLAM,PHI,
     .   FTHERM,TFAC,PLA_NL,DPDT_NL
      my_real, DIMENSION(NEL,3) :: XVEC_ETA
      my_real, DIMENSION(NEL,4) :: XVEC
      INTEGER, DIMENSION(NEL,3) :: IPOS_ETA
      INTEGER, DIMENSION(NEL,2) :: IPOS
C-----------------------------------------------
      ! VARTMP(1)   latest position of PLAS in TAB_YLD function 
      ! VARTMP(2)   latest position of PLAS in TAB_TEMP function 
      ! VARTMP(3)   latest position of TEMP in TAB_TEMP function 
      ! VARTMP(4)   latest position of TEMP in TAB_ETA function 
      ! VARTMP(5)   latest position of PLAS in TAB_ETA function 
C=======================================================================
c     Material model parameters
c-----------------------------------------------
      YOUNG   =  UPARAM(1)       ! Young modulus
      NU      =  UPARAM(2)       ! Poisson ratio
      ETA     =  UPARAM(3)       ! Thermal work coefficient
      TREF    =  UPARAM(4)       ! Reference temperature
      TINI    =  UPARAM(5)       ! Initial tempareture
      ISMOOTH =  NINT(UPARAM(6)) ! function interpolation flag
      XRATE   =  UPARAM(7)       ! strain rate abscissa factor for eta function
      XSCALE  =  UPARAM(8)       ! strain rate abscissa factor for yld function
      YSCALE  =  UPARAM(9)       ! Yld function scale factor
      
      G       =  UPARAM(11)      ! Shear modulus
      G2      =  UPARAM(12)      ! Shear modulus * 2
      G3      =  UPARAM(13)      ! Shear modulus * 3
      BULK    =  UPARAM(14)      ! Bulk modulus
      LAME    =  UPARAM(15)      ! Lame parameter
      CP      =  UPARAM(20)      ! Specific heat coefficient
c
      FUNC_YLD  = ITABLE(1)
      FUNC_TEMP = ITABLE(2)
      FUNC_ETA  = ITABLE(3)
      NDIM_YLD  = TABLE(FUNC_YLD)%NDIM
      IF (FUNC_TEMP > 0) THEN
        NDIM_TEMP  = TABLE(FUNC_TEMP)%NDIM
      ENDIF
      IF (FUNC_ETA  > 0) THEN
        NDIM_ETA   = TABLE(FUNC_ETA)%NDIM
      ENDIF
c 
      ! Maximal number of Newton iterations
      NITER = 3
c
c     Initializations
c
      DPLA(1:NEL)  = ZERO          ! Initialization of the plastic strain increment
      ET(1:NEL)    = ONE            ! Initialization of tangent stiffness factor
      DTINV  = ONE / MAX(EM20, TIMESTEP)
      DLAMIN = EM15
      ALPHA  = 0.025
      ALPHI  = ONE-ALPHA
      IF (INLOC > 0) THEN 
        DO I = 1,NEL
          UVAR(I,1)  = UVAR(I,1) + MAX(DPLANL(I),ZERO)
          PLA_NL(I)  = UVAR(I,1)
          DPDT_NL(I) = MAX(DPLANL(I),ZERO)*DTINV
        ENDDO
      ENDIF
c      
c---  self heating factor and temperatiure initialization
      IF (JTHE == 0 .and. TIME == ZERO) TEMP(1:NEL) = TINI
!
      IF (ETA > ZERO) THEN
        IF (FUNC_ETA > 0) THEN    ! scale factor function for Taylor-Quinney coefficient
          IF (INLOC == 0) THEN 
            XVEC_ETA(1:NEL,1) = EPSP(1:NEL) * XRATE 
          ELSE
            XVEC_ETA(1:NEL,1) = DPDT_NL(1:NEL) * XRATE 
          ENDIF
          IPOS_ETA(1:NEL,1) = 1
          IF (NDIM_ETA > 1) THEN
            XVEC_ETA(1:NEL,2) = TEMP(1:NEL)
            IPOS_ETA(1:NEL,2) = VARTMP(1:NEL,4)
          END IF
          IF (NDIM_ETA > 2) THEN
            IF (INLOC == 0) THEN 
              XVEC_ETA(1:NEL,3) = PLA(1:NEL)
            ELSE
              XVEC_ETA(1:NEL,3) = PLA_NL(1:NEL)
            ENDIF
            IPOS_ETA(1:NEL,3) = VARTMP(1:NEL,5)
          END IF

          CALL TABLE_VINTERP(TABLE(FUNC_ETA),NEL,NEL,IPOS_ETA,XVEC_ETA,FACT_ETA,DXDYV)          
    
          IF (NDIM_ETA > 1) VARTMP(1:NEL,4) = IPOS_ETA(1:NEL,2)
          IF (NDIM_ETA > 2) VARTMP(1:NEL,5) = IPOS_ETA(1:NEL,3)
          DO I=1,NEL
            FTHERM(I) = MIN(ETA*FACT_ETA(I), ONE)
          END DO
        ELSE
          FTHERM(1:NEL) = MIN(ETA, ONE)
        END IF 
      END IF 
c-----------------------------------------------
c     Trial stress
c-----------------------------------------------
      DO I=1,NEL
        LDAV = (DEPSXX(I) + DEPSYY(I) + DEPSZZ(I)) * LAME
        SIGNXX(I) = SIGOXX(I) + DEPSXX(I)*G2 + LDAV
        SIGNYY(I) = SIGOYY(I) + DEPSYY(I)*G2 + LDAV
        SIGNZZ(I) = SIGOZZ(I) + DEPSZZ(I)*G2 + LDAV
        SIGNXY(I) = SIGOXY(I) + DEPSXY(I)*G
        SIGNYZ(I) = SIGOYZ(I) + DEPSYZ(I)*G
        SIGNZX(I) = SIGOZX(I) + DEPSZX(I)*G
        SIGM(I)   = (SIGNXX(I) + SIGNYY(I) + SIGNZZ(I)) * THIRD
        ! deviatoric trial stress tensor
        STXX(I) = SIGNXX(I) - SIGM(I)
        STYY(I) = SIGNYY(I) - SIGM(I)
        STZZ(I) = SIGNZZ(I) - SIGM(I)
        STXY(I) = SIGNXY(I)
        STYZ(I) = SIGNYZ(I)
        STZX(I) = SIGNZX(I)
      ENDDO
c----------------------------------------------------      
c     Computation of the initial yield stress
c----------------------------------------------------      
      XVEC(1:NEL,1) = PLA(1:NEL)
      XVEC(1:NEL,2) = EPSP(1:NEL) * XSCALE
      IPOS(1:NEL,1) = VARTMP(1:NEL,1)
      IPOS(1:NEL,2) = 1
c
      CALL TABLE2D_VINTERP_LOG(TABLE(FUNC_YLD),ISMOOTH,NEL,NEL,IPOS,XVEC  ,YLD  ,HARDP,HARDR)               
c
      YLD(1:NEL)   = YLD(1:NEL)   * YSCALE
      HARDP(1:NEL) = HARDP(1:NEL) * YSCALE
      VARTMP(1:NEL,1) = IPOS(1:NEL,1)
c----------------------------------------------------      
c     Computation of temperature dependent yield stress factor from quasistatic curves
c----------------------------------------------------      
      IF (FUNC_TEMP > 0) THEN
c        XVEC(1:NEL,1) = PLA(1:NEL)
        XVEC(1:NEL,2) = TREF
        IPOS(1:NEL,1) = VARTMP(1:NEL,2)
        IPOS(1:NEL,2) = VARTMP(1:NEL,3)
        CALL TABLE_VINTERP(TABLE(FUNC_TEMP),NEL,NEL,IPOS,XVEC,YLD_TREF,DYDX)  
        VARTMP(1:NEL,2) = IPOS(1:NEL,1)     
        VARTMP(1:NEL,3) = IPOS(1:NEL,2)     
c     
        XVEC(1:NEL,2) = TEMP(1:NEL)
        CALL TABLE_VINTERP(TABLE(FUNC_TEMP),NEL,NEL,IPOS,XVEC,YLD_TEMP,DYDX)          
c     
        TFAC(1:NEL)  = YLD_TEMP(1:NEL) / YLD_TREF(1:NEL)      
        YLD(1:NEL)   = YLD(1:NEL)   * TFAC(1:NEL)      
        HARDP(1:NEL) = HARDP(1:NEL) * TFAC(1:NEL) 
      ELSE
        TFAC(1:NEL) = ONE
      END IF
c-----------------------------------------------
c     Check plasticity 
c-----------------------------------------------
      NINDX = 0
      DO I=1,NEL
        J2 = STXY(I)*STXY(I) + STYZ(I)*STYZ(I) + STZX(I)*STZX(I) 
     .     - STXX(I)*STYY(I) - STYY(I)*STZZ(I) - STZZ(I)*STXX(I)
        Q2 = THREE*J2
c
        IF (Q2 > YLD(I)**2  .and. OFF(I) == ONE) THEN
          NINDX = NINDX + 1   ! Number of the elements with plastic behaviour
          INDEX(NINDX)  = I
          SVMT(I) = SQRT(Q2)
          PHI(I)  = SVMT(I) - YLD(I)
        ENDIF
      ENDDO
c      
      !====================================================================
      ! - PLASTIC CORRECTION Newton
      !====================================================================
c
      IF (NINDX > 0) THEN
        DO II = 1, NINDX   ! Number of the element with plastic behaviour     
          I = INDEX(II)
          ! initialize incremental update
          SXX(I) = STXX(I)
          SYY(I) = STYY(I)
          SZZ(I) = STZZ(I)
          SXY(I) = STXY(I)
          SYZ(I) = STYZ(I)
          SZX(I) = STZX(I)
          SVM(I) = SVMT(I)
        END DO
c
        DO ITER = 1,NITER 
c
          DO II = 1, NINDX      
            I = INDEX(II)
c
            PHI(I)    = SVM(I) - YLD(I)
            DPHI_DLAM = G3 + HARDP(I)
            DLAM(I) = PHI(I) / DPHI_DLAM
            DPLA(I) = DPLA(I) + DLAM(I)
          END DO
          
          ! Update Yld and Hardp with new plastic strain and strain rate
          
          DO II = 1, NINDX 
            I = INDEX(II)
            XVEC(II,1) = PLA(I) + DPLA(I)
            XVEC(II,2) = EPSP(I)
            IPOS(II,1) = VARTMP(I,1)
            IPOS(II,2) = 1
          ENDDO
c
          CALL TABLE2D_VINTERP_LOG(TABLE(FUNC_YLD),ISMOOTH,NEL,NINDX,IPOS,XVEC,YLD_I,HARDP_I,HARDR_I)
c
          DO II = 1, NINDX 
            I = INDEX(II)
            VARTMP(I,1) = IPOS(II,1)
            HARDP(I) = HARDP_I(II)  *YSCALE*TFAC(I)
            YLD(I)   = MIN(YLD_I(II)*YSCALE*TFAC(I), SVMT(I))
          ENDDO
c                    
          ! deviatoric stress update
c
          DO II = 1, NINDX 
            I = INDEX(II)
c
            R = G3 * DPLA(I) / SVMT(I)
            SXX(I) = STXX(I) - R * SXX(I)
            SYY(I) = STYY(I) - R * SYY(I)
            SZZ(I) = STZZ(I) - R * SZZ(I)
            SXY(I) = STXY(I) - R * SXY(I)
            SYZ(I) = STYZ(I) - R * SYZ(I)
            SZX(I) = STZX(I) - R * SZX(I)

            J2 = SXY(I)*SXY(I) + SYZ(I)*SYZ(I) + SZX(I)*SZX(I) 
     .         - SXX(I)*SYY(I) - SYY(I)*SZZ(I) - SZZ(I)*SXX(I)
            SVM(I) = SQRT(THREE*J2)
            PHI(I) = SVM(I) - YLD(I)
         ENDDO              
c            
        END DO   !  ITER = 1,NITER
c        
        DO II = 1, NINDX 
          I = INDEX(II)
          SIGNXX(I) = SXX(I) + SIGM(I)
          SIGNYY(I) = SYY(I) + SIGM(I)
          SIGNZZ(I) = SZZ(I) + SIGM(I)
          SIGNXY(I) = SXY(I)
          SIGNYZ(I) = SYZ(I)
          SIGNZX(I) = SZX(I)
c
          PLA(I) = PLA(I) + MAX(DPLA(I), ZERO)
          ET(I)  = HARDP(I) / (HARDP(I) + YOUNG)
        ENDDO                      
c
        ! Update the temperature or heat energy due to plastic work
        IF (ETA > ZERO .AND. INLOC == 0) THEN 
          IF (JTHE /= 0 .AND. JLAG /= 0) THEN
            DO II = 1, NINDX 
              I = INDEX(II)
              FHEAT(I) = FHEAT(I) + FTHERM(I)*YLD(I)*DPLA(I) * VOLUME(I)
            ENDDO
          ELSE IF (CP > ZERO) THEN
            DO II = 1, NINDX 
              I = INDEX(II)
              TEMP(I) = TEMP(I) + FTHERM(I)*YLD(I)*DPLA(I) / (CP*RHO(I))
            ENDDO
          ENDIF
        ENDIF
c
      END IF     !  NINDX > 0  (plasticity)
!-----------------------------------------------------------------------------
      ! if EOS is used, material law calculates only deviatoric stress tensor
      !                 sound speed depends on pressure derivative over volume change
      !                 calculated in EOS
      IF (IEOS > 0) THEN
        SIGNXX(1:NEL) = SIGNXX(1:NEL) - SIGM(1:NEL)
        SIGNYY(1:NEL) = SIGNYY(1:NEL) - SIGM(1:NEL)
        SIGNZZ(1:NEL) = SIGNZZ(1:NEL) - SIGM(1:NEL)
        SIGNXY(1:NEL) = SIGNXY(1:NEL)
        SIGNYZ(1:NEL) = SIGNYZ(1:NEL)
        SIGNZX(1:NEL) = SIGNZX(1:NEL)
        SOUNDSP(1:NEL)= SQRT((DPDM(1:NEL) + FOUR_OVER_3*G) / RHO(1:NEL))
      ELSE 
        SOUNDSP(1:NEL)= SQRT((BULK + FOUR_OVER_3*G) / RHO(1:NEL))
      END IF                
c
      SIGY(1:NEL) = YLD(1:NEL)
      EPSP(1:NEL) = ALPHA*DPLA(1:NEL)*DTINV + ALPHI*EPSP(1:NEL)
      IF (INLOC > 0 .AND. ETA > ZERO) THEN 
        IF (JTHE /= 0 .AND. JLAG /= 0) THEN
          DO II = 1, NEL
            FHEAT(I) = FHEAT(I) + FTHERM(I)*YLD(I)*MAX(DPLANL(I),ZERO) * VOLUME(I)
          ENDDO
        ELSE IF (CP > ZERO) THEN
          DO I = 1,NEL     
            TEMP(I) = TEMP(I) + FTHERM(I)*YLD(I)*MAX(DPLANL(I),ZERO) / (CP*RHO(I))
          ENDDO
        ENDIF
      ENDIF
c-----------
      RETURN
      END SUBROUTINE SIGEPS109
