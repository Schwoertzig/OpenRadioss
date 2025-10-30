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
      module hm_read_mat251_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief General Multimaterial Law 251
!! \details if needed, more details can be added here      
      subroutine hm_read_mat251( &
                   npropmi,  npropm , ipm     ,pm     ,mat_id  ,titr  , &
                   lsubmodel,mtag   ,matparam ,unitab)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use elbuftag_mod
      use unitab_mod , only : unit_type_
      use constant_mod , only : em01
      use message_mod
      use submodel_mod , only : submodel_data, nsubmod
      use matparam_def_mod , only : matparam_struct_
      use names_and_titles_mod , only : nchartitle
      use multimat_param_mod , only : m20_discrete_fill
      use precision_mod , only : wp
      use constant_mod , only : one
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "units_c.inc"
    COMMON /COM01/ NMULT
    INTEGER NMULT
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      INTEGER,INTENT(IN) :: NPROPM, NPROPMI
      INTEGER, INTENT(INOUT) :: IPM(NPROPMI)
      INTEGER, INTENT(IN) :: MAT_ID
      real(kind=WP), INTENT(INOUT) :: PM(NPROPM)
      TYPE(MLAW_TAG_), INTENT(INOUT) :: MTAG
      CHARACTER(LEN=NCHARTITLE),INTENT(IN) ::  TITR
      TYPE(SUBMODEL_DATA),INTENT(IN) :: LSUBMODEL(NSUBMOD)
      TYPE(MATPARAM_STRUCT_) ,INTENT(INOUT) :: MATPARAM
      TYPE (UNIT_TYPE_),INTENT(IN) :: UNITAB
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      integer :: imat,nmat
      integer :: submat_id
      real(kind=wp) :: frac_vol
      logical :: is_encrypted, is_available
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      is_encrypted = .false.
      is_available = .false.
      call hm_option_is_encrypted(is_encrypted)

      call hm_get_intv('NIP',nmat,is_available,lsubmodel)
      matparam%multimat%nb=nmat
      if(.not.allocated(matparam%multimat%vfrac))allocate(matparam%multimat%vfrac(nmat))
      if(.not.allocated(matparam%multimat%mid))  allocate(matparam%multimat%mid(nmat))

      do imat=1,nmat
        call hm_get_int_array_index('MAT_ID_ARRAY',submat_id,imat,is_available,lsubmodel)
        call hm_get_float_array_index('VOL_FRAC_ARRAY',frac_vol,imat,is_available,lsubmodel,unitab)
        matparam%multimat%vfrac(imat)=frac_vol
        matparam%multimat%mid(imat)=submat_id
      enddo
      
      do imat=1,nmat
         ipm(20 + imat) = matparam%multimat%mid(imat)
         pm(20 + imat) =matparam%multimat%vfrac(imat)
      enddo

      pm(1) = one
      pm(89) =one

      pm(20)=nmat+em01
      nmult=max(nmult,nmat)

      ! output
      write(iout, 900) trim(titr), mat_id, 251
      write(iout,1000) nmat
      if(is_encrypted)then
        write(iout,'(5x,a,//)')'confidential data'
      else
        do imat=1,nmat
          write(iout,1010) matparam%multimat%mid(imat),matparam%multimat%vfrac(imat)
        enddo
      endif
      
      !definition of internal variables (elementary storage)
      mtag%l_frac = 1   ! percentage of phase
      mtag%g_rho = 1
      mtag%l_rho = 1


      ! MATPARAM keywords
      CALL INIT_MAT_KEYWORD(MATPARAM,"INCOMPRESSIBLE")

      ! Material compatibility with /EOS option
      CALL INIT_MAT_KEYWORD(MATPARAM,"EOS")

      ! EOS/Thermo keyword for pressure treatment in elements
      CALL INIT_MAT_KEYWORD(MATPARAM,"HYDRO_EOS")

      ! Properties compatibility
      CALL INIT_MAT_KEYWORD(MATPARAM,"SOLID_ISOTROPIC")

      return
!--------------------------------
  900 format(/&
      5x,a,/,&
      5x,'material number. . . . . . . . . . . . . . .=',i10/,&
      5x,'material law . . . . . . . . . . . . . . . .=',i10/)
 1000 format(&
      5x,'  multimat law /mat/mmale   ',/,&
      5x,'  -----------------------  ',/&
      5x,'  number of materials. . . . . . .=',i10//)
 1010 format(&
     & 5x,'  material id', i10, ' ; volume fraction',1pg20.13/)
!--------------------------------
      return
      end subroutine hm_read_mat251
      end module hm_read_mat251_mod
