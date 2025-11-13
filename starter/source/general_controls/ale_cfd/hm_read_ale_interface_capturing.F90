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
      module HM_READ_ALE_INTERFACE_CAPTURING_MOD
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Reader for keyword /ALE/INTERFACE_CAPTURING/...
!! \details
        subroutine HM_READ_ALE_INTERFACE_CAPTURING(LSUBMODEL, UNITAB)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
         USE MESSAGE_MOD
         USE HM_OPTION_READ_MOD
         USE NAMES_AND_TITLES_MOD , ONLY : NCHARFIELD
         USE ALE_MOD , ONLY : ALE
         USE SUBMODEL_MOD , ONLY : NSUBMOD, SUBMODEL_DATA
         USE UNITAB_MOD , ONLY : UNIT_TYPE_
         USE CONSTANT_MOD , only : zero,two,one
         USE PRECISION_MOD , ONLY : WP
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      TYPE(SUBMODEL_DATA), DIMENSION(NSUBMOD), INTENT(IN) :: LSUBMODEL
      TYPE(UNIT_TYPE_), INTENT(IN) :: UNITAB
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: NALE_INTER_CAPTURING !< number of defined option
      CHARACTER(LEN=NCHARFIELD) :: KEY
      LOGICAL :: IS_AVAILABLE
      REAL(kind=WP) :: BETA,KAPPA
      INTEGER :: II
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Precondition
! ---------------------------------------------------------------------------------------------------------------------
! None
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      ALE%INTER_CAPTURING%IS_DEFINED = 0
      ALE%INTER_CAPTURING%IMETH = 2

      CALL HM_OPTION_COUNT('/ALE/INTERFACE_CAPTURING', NALE_INTER_CAPTURING)
      IF(NALE_INTER_CAPTURING > 1)THEN
        !Warning message : multiple definition has no sens. Last read definition will be used
        CALL ANCMSG(MSGID = 140, MSGTYPE=MSGWARNING, ANMODE=ANINFO)
      ENDIF

      CALL HM_OPTION_START('/ALE/INTERFACE_CAPTURING')

      DO II = 1, NALE_INTER_CAPTURING

         CALL HM_OPTION_READ_KEY(LSUBMODEL, KEYWORD3 = KEY)
         BETA = ZERO
         KAPPA = ZERO

         !------------------------------------------------------!
         !     /ALE/INTERFACE_CAPTURING/UPWIND                  !
         !------------------------------------------------------!
         IF (KEY(1:6) == 'UPWIND') THEN
            ALE%INTER_CAPTURING%IS_DEFINED = 1
            ALE%INTER_CAPTURING%IMETH = 1

         !------------------------------------------------------!
         !     /ALE/INTERFACE_CAPTURING/MUSCL                   !
         !------------------------------------------------------!
         ELSEIF (KEY(1:5) == 'MUSCL') THEN
            ALE%INTER_CAPTURING%IS_DEFINED = 1
            ALE%INTER_CAPTURING%IMETH = 2
            CALL HM_GET_FLOATV('ALE_MUSCL_BETA', BETA, IS_AVAILABLE, LSUBMODEL, UNITAB)
            IF(BETA < ZERO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO,C1='MUSCL',C2='BETA',C3='PARAMETER MUST BE POSITIVE')
             END IF
            IF(BETA > TWO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO,C1='MUSCL',C2='BETA',C3='PARAMETER SHOULD BE LOWER THEN 2.0')
             END IF
             IF(BETA == ZERO) BETA = TWO

         !------------------------------------------------------!
         !     /ALE/INTERFACE_CAPTURING/SUPERBEE                !
         !------------------------------------------------------!
         ELSEIF (KEY(1:8) == 'SUPERBEE' .OR. KEY(1:9) == 'SUPER-BEE') THEN
            ALE%INTER_CAPTURING%IS_DEFINED = 1
            ALE%INTER_CAPTURING%IMETH = 3
            CALL HM_GET_FLOATV('ALE_SUPERBEE_BETA', BETA, IS_AVAILABLE, LSUBMODEL, UNITAB)
            IF(BETA < ZERO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO,C1='SUPERBEE',C2='BETA',C3='PARAMETER MUST BE POSITIVE')
             END IF
            IF(BETA > TWO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO, &
               C1='SUPERBEE',C2='BETA',C3='PARAMETER SHOULD BE LOWER THEN 2.0')
             END IF
             IF(BETA == ZERO) BETA = TWO

         !------------------------------------------------------!
         !     /ALE/INTERFACE_CAPTURING/HYPER-C                 !
         !------------------------------------------------------!
         ELSEIF (KEY(1:7) == 'HYPER-C' .OR. KEY(1:6) == 'HYPERC' .OR. KEY(1:7) == 'HYPER_C') THEN
           ALE%INTER_CAPTURING%IS_DEFINED = 1
           ALE%INTER_CAPTURING%IMETH = 4
           CALL HM_GET_FLOATV('ALE_HYPERC_BETA', BETA, IS_AVAILABLE, LSUBMODEL, UNITAB)
           CALL HM_GET_FLOATV('ALE_HYPERC_KAPPA', KAPPA, IS_AVAILABLE, LSUBMODEL, UNITAB)
           IF(BETA < ZERO)THEN
              CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO,C1='SUPERBEE',C2='BETA',C3='PARAMETER MUST BE POSITIVE')
            END IF
           IF(BETA > TWO)THEN
              CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO, &
              C1='HYPER-C',C2='BETA',C3='PARAMETER SHOULD BE LOWER THEN 2.0')
            END IF
            IF(BETA == ZERO) BETA = TWO
            IF(KAPPA == ZERO) KAPPA = ONE

          !------------------------------------------------------!
          !     /ALE/INTERFACE_CAPTURING/DEPRES-LAGOUTIERE       !
          !------------------------------------------------------!
         ELSEIF (KEY(1:17) == 'DEPRES-LAGOUTIERE' .OR. KEY(1:2) == 'DL') THEN
            ALE%INTER_CAPTURING%IS_DEFINED = 1
            ALE%INTER_CAPTURING%IMETH = 5
            CALL HM_GET_FLOATV('ALE_DL_BETA', BETA, IS_AVAILABLE, LSUBMODEL, UNITAB)
            CALL HM_GET_FLOATV('ALE_DL_KAPPA', KAPPA, IS_AVAILABLE, LSUBMODEL, UNITAB)
            IF(BETA < ZERO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO,C1='SUPERBEE',C2='BETA',C3='PARAMETER MUST BE POSITIVE')
             END IF
            IF(BETA > TWO)THEN
               CALL ANCMSG(MSGID=141,MSGTYPE=MSGWARNING,ANMODE=ANINFO, &
               C1='DEPRES-LAGOUTIERE',C2='BETA',C3='PARAMETER SHOULD BE LOWER THEN 2.0')
             END IF
             IF(BETA == ZERO) BETA = TWO
             IF(KAPPA == ZERO) KAPPA = ONE

         ELSE
           !invalid keyword : error message managed by Reader itself
     
         ENDIF

         ALE%INTER_CAPTURING%BETA = BETA
         ALE%INTER_CAPTURING%KAPPA = KAPPA

      ENDDO

       ! POST-VERIFICATION
       IF(ALE%INTER_CAPTURING%IS_DEFINED /= 0 .AND. ALE%INTER_CAPTURING%IS_DEFINED /= 1) ALE%INTER_CAPTURING%IS_DEFINED = 0
       IF(ALE%INTER_CAPTURING%IMETH <= 0 .OR. ALE%INTER_CAPTURING%IMETH >= 6) ALE%INTER_CAPTURING%IMETH = 2

      END SUBROUTINE HM_READ_ALE_INTERFACE_CAPTURING
      END MODULE HM_READ_ALE_INTERFACE_CAPTURING_MOD
