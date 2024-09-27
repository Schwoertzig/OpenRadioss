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
      SUBROUTINE HM_READ_EBCS_PROPERGOL(IGRSURF, MULTI_FVM, UNITAB, ID, TITR, UID, LSUBMODEL,  NSURF, NUMNOD, EBCS)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use ebcs_mod
      use unitab_mod
      use message_mod
      use multi_fvm_mod
      use groupdef_mod
      use submodel_mod
      use names_and_titles_mod , only : nchartitle, ncharkey
      use constant_mod , only : zero, one
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
#include      "units_c.inc"
#include      "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer,intent(in) :: nsurf, numnod !< array sizes
      type (unit_type_),intent(in) ::unitab
      integer id,uid
      type (multi_fvm_struct), intent(inout) :: multi_fvm
      type (surf_)   ,target,  dimension(nsurf)   :: igrsurf
      character(len=nchartitle), intent(in) :: titr
      type(submodel_data) lsubmodel(nsubmod)
      logical is_available,is_encrypted
      type(t_ebcs_propergol), intent(inout) :: ebcs
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      integer i,isu,sens,monvol,surf,ipres,irho,nod,iad,j,k1,k2,nseg,iener,ivx,ivy,ivz,ialpha
      integer liste(numnod), imat,ivel_typ,u_ialpha,u_irho,u_ipres,iflagunit,flag_fmt,check_cumul_vf(2)
      my_real c,pres,rho,lcar,r1,r2,ener,vx,vy,vz, alpha, param_a, param_n, param_q, param_rho0s
      character mess*40,chain*9, chain1*64
      character(len=ncharkey) :: key
      logical found
      integer, dimension(:), pointer :: ingr2usr
      integer, external :: ngr2usr
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      
      iflagunit=0
      do j=1,unitab%nunits
        if (unitab%unit_id(j) == uid) then
          iflagunit = 1
          exit
        endif
      enddo
      if (uid/=0.and.iflagunit==0) then
        call ancmsg(msgid=659,anmode=aninfo,msgtype=msgerror,i2=uid,i1=id,c1='ebcs',c2='ebcs',c3=titr)
      endif
    
      call hm_option_is_encrypted(is_encrypted)
      call hm_get_intv('entityid', surf, is_available,lsubmodel)
      call hm_get_floatv('param_a', param_a, is_available,lsubmodel,unitab)
      call hm_get_floatv('param_n', param_n, is_available,lsubmodel,unitab)
      call hm_get_floatv('param_q', param_q, is_available,lsubmodel,unitab)
      call hm_get_floatv('rho0s', param_rho0s, is_available,lsubmodel,unitab)

      if(param_a == zero)param_a = one
      if(param_n < zero)param_n = zero
      if(param_q < zero)param_q = zero
      if(param_rho0s == zero)param_rho0s = one

      ebcs%a = param_a
      ebcs%n = param_n
      ebcs%q = param_q
      ebcs%rho0s = param_rho0s
      ebcs%has_ielem = .true.

      if(multi_fvm%is_used)then
        ebcs%is_multifluid = .true.
      endif

       ebcs%fvm_inlet_data%func_vel(1:3) = -1
       ebcs%fvm_inlet_data%val_vel(1:3) = zero
       ebcs%fvm_inlet_data%formulation = 2
       ebcs%fvm_inlet_data%vector_velocity = 1
       do imat = 1, 21 !  multi_fvm%nbmat       -> init from nbmat+1,21 to avoid uninit values transmitted to starter
          ebcs%fvm_inlet_data%func_alpha(imat) = -1
          ebcs%fvm_inlet_data%func_rho(imat)   = -1
          ebcs%fvm_inlet_data%func_pres(imat)  = -1
          ebcs%fvm_inlet_data%val_alpha(imat)  = zero
          ebcs%fvm_inlet_data%val_rho(imat)    = zero
          ebcs%fvm_inlet_data%val_pres(imat)   = zero
       enddo

         isu=0
         ingr2usr => igrsurf(1:nsurf)%id
         if (surf /= 0) isu=ngr2usr(surf,ingr2usr,nsurf)
         nseg=0
         if (isu /= 0) nseg=igrsurf(isu)%nseg
         if(surf == 0)then
            ierr=ierr+1
            write(istdo,'(6X,A)')' ** A SURFACE SHOULD BE INPUT'
            write(iout, '(6X,A)')' ** A SURFACE SHOULD BE INPUT'
         elseif(isu == 0)then
            ierr=ierr+1
            write(istdo,*)' ** ERROR SURFACE NOT FOUND, ID=',SURF
            write(iout,*) ' ** ERROR SURFACE NOT FOUND, ID=',SURF
         elseif(nseg == 0)then
            ierr=ierr+1
            write(istdo,*)' ** ERROR EMPTY SURFACE',SURF
            write(iout,*) ' ** ERROR EMPTY SURFACE',SURF
         endif

         ebcs%nb_elem = nseg

         write(iout,1001)id, trim(titr)
         write(iout,1118)surf,nseg,param_rho0s, param_q
         write(iout,1201)param_a,param_n
         
      
!-----------
      return
!-----------
 1001 format( //'PROPERGOL EBCS NUMBER. . . . :',I8,1X,A)

 1118 format(&
              '    ON SURFACE  . . . . . . . . . . . . . . . ',I8,/,&
              '    NUMBER OF SEGMENTS FOUND. . . . . . . . . ',I8,/,&
              '    PROPERGOL DENSITY . . . . . . . . . . . . ',E20.12,/,&
              '    PROPERGOL HEAT OF COMBUSTION. . . . . . . ',E20.12)
 1201 format( '      --- COMBUSION MODEL : VIEILLE''S LAW'    ,/,&
              '      VIEILLE PARAMETER A . . . . . . . . . . ',E20.12,/,&
              '      VIEILLE PARAMETER N . . . . . . . . . . ',E20.12 )

      end subroutine hm_read_ebcs_propergol

