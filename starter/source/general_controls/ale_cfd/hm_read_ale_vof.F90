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
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
      SUBROUTINE HM_READ_ALE_VOF(LSUBMODEL,N2D)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      USE HM_OPTION_READ_MOD
      USE SUBMODEL_MOD
      USE MESSAGE_MOD
      USE ALE_MOD , ONLY : ALE
       USE ALEMUSCL_MOD , only : ALEMUSCL_Param
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
       implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      TYPE(SUBMODEL_DATA), DIMENSION(NSUBMOD), INTENT(IN) :: LSUBMODEL
      INTEGER, INTENT(IN) :: N2D
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER :: NALEVOF, NALEMUSCL
      INTEGER :: IFORM
      LOGICAL :: IS_AVAILABLE
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      CALL HM_OPTION_COUNT('/ALE/VOF', NALEVOF)


      IF (NALEVOF > 0) THEN
         CALL HM_OPTION_START('/ALE/VOF')
         CALL HM_OPTION_NEXT()
         CALL HM_GET_INTV('IFORM', IFORM, IS_AVAILABLE, LSUBMODEL)
         IF (IFORM < 0 .OR. IFORM > 1)THEN
           IFORM = 0
           CALL ANCMSG(MSGID   = 1139, MSGTYPE = MSGERROR, ANMODE = ANINFO, C1 = 'IFORM', I1 = IFORM )
         ENDIF
         IF(IFORM == 1)THEN
           IF(N2D == 0)THEN
               CALL ANCMSG(MSGID   = 1140, MSGTYPE = MSGERROR, ANMODE = ANINFO,  &
               C1 = 'NOT COMPATIBLE WITH 3D ANALYSIS', I1 = IFORM )
               IFORM = 0
           END IF
           ALE%VOF%IS_DEFINED = IFORM
           IF (IFORM == 1 .AND. ALEMUSCL_Param%IALEMUSCL /= 0)THEN
               ALEMUSCL_Param%IALEMUSCL = 0  ! LAW51 interface reconstruction :either UPWIND, or MUSCL or VOF. No possible combination
               CALL HM_OPTION_COUNT('/ALE/MUSCL', NALEMUSCL)
               IF (NALEMUSCL > 0) THEN
                  CALL ANCMSG(MSGID   = 1140, MSGTYPE = MSGWARNING, ANMODE = ANINFO,  &
                  C1 = 'MUSCL REPLACED WITH VOF (LAW51 INTERFACE RECONSTRUCTION)', I1 = IFORM )
               ENDIF
               IFORM = 0
           END IF
         END IF
      ENDIF




      END 
