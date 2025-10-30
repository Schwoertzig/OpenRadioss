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
      module ale_allocations_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
        subroutine ale_allocations(numnod,numelq,numels,numeltg,glob_therm)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
!  [ the module names in use must be in uppercase for now, it will change later]
!  [ ONLY is mandatory, note the space before the ,]
          use ale_mod , only : ale
          use glob_therm_mod , only : glob_therm_
          use message_mod
          use restmod , only : dflow,vflow,wflow ! foam material law
          use ALEFVM_MOD , only : ALEFVM_Buffer, ALEFVM_Param
          use constant_mod , only : zero
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer,intent(in) :: numnod                      !< number of nodes
          type(glob_therm_) ,intent(inout) :: glob_therm    !< global data structure for heat equation
          integer,intent(in) :: numels,numelq,numeltg       !< number of solid elements
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: iale,ieuler,ialelag
          integer :: nmult
          integer :: is_defined_law20
          integer :: stat
! ----------------------------------------------------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
! [ external functions must be kept to minimum ]
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      IALE = ALE%GLOBAL%IALE
      IEULER = ALE%GLOBAL%IEULER
      IALELAG = ALE%GLOBAL%IALELAG
      NMULT = ALE%MULTIMAT%NMULT
      IS_DEFINED_LAW20 = ALE%MULTIMAT%IS_DEFINED_LAW20
      ALE%GLOBAL%SNALE  = MAX(IALE,IEULER,IALELAG)*NUMNOD !nodal tags (ALE,EULER,LAG, ...)
      ALE%GLOBAL%SIELVS = 6*NUMELS+MAX(IALE,GLOB_THERM%ITHERM,IEULER,IALELAG)* (4 * NUMELQ + 3 * NUMELTG) ! adjacent elems (soze)
      ALE%MULTIMAT%VOF%SIFILL = NMULT*NUMNOD*IS_DEFINED_LAW20 !VOF law20
      ALE%MULTIMAT%VOF%SIMS   = NMULT*NUMNOD*IS_DEFINED_LAW20 ! VOF law20
      ALLOCATE(ALE%MULTIMAT%VOF%IFILL(ALE%MULTIMAT%VOF%SIFILL),STAT=stat)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268, ANMODE=ANINFO, MSGTYPE=MSGERROR, C1='IFILL')
      ALLOCATE(ALE%MULTIMAT%VOF%IMS(ALE%MULTIMAT%VOF%SIMS),STAT=stat)
      IF(STAT /= 0) CALL ANCMSG(MSGID=268, ANMODE=ANINFO, MSGTYPE=MSGERROR, C1='IMS')
      IF(ALE%MULTIMAT%VOF%SIFILL > 0) ALE%MULTIMAT%VOF%IFILL(:) = 0
      IF(ALE%MULTIMAT%VOF%SIMS   > 0) ALE%MULTIMAT%VOF%IMS(:) = 0

      ALLOCATE(DFLOW(3*NUMNOD*IALELAG)      ,STAT=stat)
      ALLOCATE(VFLOW(3*NUMNOD*IALELAG)      ,STAT=stat)
      ALLOCATE(WFLOW(3*NUMNOD*IALELAG)      ,STAT=stat)

      IF(IALELAG > 0) THEN
        DFLOW(:) = ZERO
        VFLOW(:) = ZERO
        WFLOW(:) = ZERO
      ENDIF

      IF(ALEFVM_Param%IEnabled > 0)THEN
        ALLOCATE(ALEFVM_Buffer%FCELL(6,NUMELS)   ,STAT=stat)
        ALEFVM_Buffer%FCELL(:,:) = ZERO
      ENDIF
      
! ----------------------------------------------------------------------------------------------------------------------
        end subroutine ale_allocations
      end module ale_allocations_mod



