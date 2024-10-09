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
      !||====================================================================
      !||    powder_burn   ../common_source/eos/powder_burn.F
      !||--- called by ------------------------------------------------------
      !||    eosmain       ../common_source/eos/eosmain.F
      !||====================================================================
      SUBROUTINE POWDER_BURN    (
     1                            IFLAG, NEL  ,PM   ,OFF  , EINT  ,MU          ,MU2 ,
     2                            ESPE , DVOL ,DF   ,VNEW , MAT   ,PSH         ,IPM ,
     3                            PNEW , DPDM ,DPDE ,THETA, ECOLD ,MAT_PARAM)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C This subroutine contains numerical solving
C of Deflagration EoS (experimental)
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE MATPARAM_DEF_MOD, ONLY : MATPARAM_STRUCT_
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "comlock.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "mvsiz_p.inc"
#include      "param_c.inc"
#include      "com04_c.inc"
#include      "com06_c.inc"
#include      "com08_c.inc"
#include      "vect01_c.inc"
#include      "scr06_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER MAT(NEL), IFLAG, NEL
      my_real,INTENT(IN) :: PM(NPROPM,NUMMAT),OFF(NEL),MU(NEL),MU2(NEL),ESPE(NEL),DVOL(NEL),DF(NEL),VNEW(NEL)
      my_real,INTENT(INOUT) :: THETA(NEL),ECOLD(NEL),PNEW(NEL),DPDM(NEL),DPDE(NEL),EINT(NEL)
      INTEGER,INTENT(IN) :: IPM(NPROPMI, NUMMAT)
      TYPE(MATPARAM_STRUCT_), INTENT(IN) :: MAT_PARAM !< material data structure
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER I, MX
      my_real :: bulk,P0,PSH,DD,EG,Gr,CC,ALPHA,FSCALE_B,FSCALE_P,FSCALE_G,FSCALE_RHO, C1,C2
      my_real :: SPH, TFEXTT
      INTEGER :: funcb,funcg
C-----------------------------------------------
C   S o u r c e  L i n e s
C-----------------------------------------------
      IF(IFLAG == 0) THEN

      BULK       = MAT_PARAM%EOS%UPARAM(01)
      P0         = MAT_PARAM%EOS%UPARAM(02)
      PSH        = MAT_PARAM%EOS%UPARAM(03)
      DD         = MAT_PARAM%EOS%UPARAM(04)
      EG         = MAT_PARAM%EOS%UPARAM(05)
      Gr         = MAT_PARAM%EOS%UPARAM(06)
      CC         = MAT_PARAM%EOS%UPARAM(07)
      ALPHA      = MAT_PARAM%EOS%UPARAM(08)
      FSCALE_B   = MAT_PARAM%EOS%UPARAM(09)
      FSCALE_P   = MAT_PARAM%EOS%UPARAM(10)
      FSCALE_G   = MAT_PARAM%EOS%UPARAM(11)
      FSCALE_RHO = MAT_PARAM%EOS%UPARAM(12)
      C1         = MAT_PARAM%EOS%UPARAM(13)
      C2         = MAT_PARAM%EOS%UPARAM(14)
      funcb      = MAT_PARAM%EOS%FUNC(1)
      funcg      = MAT_PARAM%EOS%FUNC(2)

       DO I=1,NEL
         DPDM(I) = zero   !total derivative
         DPDE(I) = zero   !partial derivative
       ENDDO

       DO I=1,NEL
         ECOLD(I)=-THREE100*SPH
       ENDDO

      ELSEIF(IFLAG == 1) THEN

       TFEXTT = ZERO
       DO I=1,NEL
         PNEW(I) = ZERO*OFF(I)
         EINT(I) = ZERO
         PNEW(I) = PNEW(I) + PSH
         TFEXTT  = TFEXTT-DVOL(I)*PSH
       ENDDO
       
#include "atomic.inc"
       TFEXT = TFEXT + TFEXTT
#include "atomend.inc"
C
C------------------------
C     COMPUTE TEMPERATURE
C------------------------

       SPH = PM(69,MAT(1))

       DO I=1,NEL
         IF(OFF(I) < ONE .OR. SPH <= ZERO) CYCLE
         THETA(I) = (DF(I)*EINT(I)/VNEW(I)-ECOLD(I))/SPH
         THETA(I) = MAX(ZERO,THETA(I))
       ENDDO
      ELSEIF (IFLAG == 2) THEN

         DO I=1, NEL
            IF (VNEW(I) > ZERO) THEN
               PNEW(I) = zero
               DPDM(I) = zero !total derivative
               DPDE(I) = zero !partial derivative
            ENDIF
         ENDDO
      ENDIF
C------------------------      
      RETURN
      END
