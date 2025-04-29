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
      !||====================================================================
      !||    hm_read_eos_compaction_tab_mod   ../starter/source/materials/eos/hm_read_eos_compaction_tab.F90
      !||--- called by ------------------------------------------------------
      !||    hm_read_eos                   ../starter/source/materials/eos/hm_read_eos.F
      !||====================================================================
      module hm_read_eos_compaction_tab_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Compaction EoS Reader (/EOS/COMPACTION)
!! \details  RHOI = PM(89)   -> provided by /MAT
!! \details  RHOR = PM(01)   -> provided by /MAT (can be erased by EOS if present : obsolete)
!! \details  PM(31) = P(MU0,E0) -> will be used to initialize diagonal of stress tensor SIG(1:3,*)
      !||====================================================================
      !||    hm_read_eos_compaction_tab   ../starter/source/materials/eos/hm_read_eos_compaction_tab.F90
      !||--- called by ------------------------------------------------------
      !||    hm_read_eos               ../starter/source/materials/eos/hm_read_eos.F
      !||--- calls      -----------------------------------------------------
      !||    ancmsg                    ../starter/source/output/message/message.F
      !||    finter                    ../starter/source/tools/curve/finter.F
      !||    hm_get_floatv             ../starter/source/devtools/hm_reader/hm_get_floatv.F
      !||    hm_get_intv               ../starter/source/devtools/hm_reader/hm_get_intv.F
      !||    hm_option_is_encrypted    ../starter/source/devtools/hm_reader/hm_option_is_encrypted.F
      !||--- uses       -----------------------------------------------------
      !||    elbuftag_mod              ../starter/share/modules1/elbuftag_mod.F
      !||    message_mod               ../starter/share/message_module/message_mod.F
      !||    submodel_mod              ../starter/share/modules1/submodel_mod.F
      !||====================================================================
      subroutine hm_read_eos_compaction_tab(iout,pm,unitab,lsubmodel,uid,eos_tag,ieos,npropm, &
                                            maxeos,eos_param,ntable,table)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use message_mod
      use unitab_mod , only : unit_type_
      use submodel_mod , only : nsubmod, submodel_data
      use elbuftag_mod , only : eos_tag_
      use constant_mod , only : zero, em20, em12, half, two_third, one, two, three, three100, ep10, ep20
      use eos_param_mod , only : eos_param_
      use table_mod , only : ttable
      use eos_table_copy_mod , only : eos_table_copy
      use table_mat_vinterp_mod , only : table_mat_vinterp
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer,intent(in) :: npropm, maxeos  !< array sizes
      type (unit_type_),intent(in) ::unitab !< data structure for units (/UNIT)
      integer, intent(in) :: iout !< file units
      my_real, intent(inout) :: pm(npropm)  !< data structure for material laws
      type(submodel_data), dimension(nsubmod), intent(in) :: lsubmodel !< data structure for sumobeling method (//SUBMODEL)
      integer,intent(in) :: uid
      type(eos_tag_),dimension(0:maxeos) ,intent(inout) :: eos_tag !< data structure for EoS
      integer,intent(in) :: ieos !< EoS (internal) identifier
      type(eos_param_), intent(inout) :: eos_param !< eos data structure (specific parameters)
      integer, intent(in) :: ntable
      type(ttable) ,dimension(ntable) ,intent(in) :: table
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
      my_real :: p0, psh, rhoi,rhor
      my_real :: ssp0
      logical :: is_encrypted, is_available
      integer :: P_FUNC_ID, C_FUNC_ID !< user function identifer
      integer :: PLASEXP
      integer :: iform
      integer :: ierror
      integer :: npt
      my_real :: PSCALE, CSCALE
      my_real :: RHO_TMD, C_SOLID
      my_real :: density_unit
      my_real :: rhomax_plastic

      my_real :: x1scale, x2scale, x3scale, x4scale
      my_real :: x2vect(2),x3vect(2),x4vect(2)
      my_real :: fscale(2)

      my_real :: slope_end

      my_real :: puser, ff, df, rho_, tol, residu  !< variable for newton iteration
      integer :: niter, iter !< variable for newton iteration

      my_real :: x1,x2,y1,y2,dydx,dx,dy
      my_real, dimension(1,1) :: xvec1 !<temporary array for table interpolation

      integer :: vartmp(1,2)
      my_real :: slope(1)
      my_real :: yy(1)

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      is_encrypted = .false.
      is_available = .false.
      iform=1
      vartmp(1,1:2) = 1

      eos_tag(ieos)%g_mu = 1
      eos_tag(ieos)%l_mu = 1
      eos_tag(ieos)%nvar = 3   !  --> elbuf%bufly%eos%nvar      (LAMBDA, C*C, P)
      eos_tag(ieos)%nvartmp = 2
           
      call hm_option_is_encrypted(is_encrypted)

      call hm_get_intv('P_FUNC', P_FUNC_ID, is_available,lsubmodel)
      call hm_get_floatv('PSCALE', Pscale, is_available,lsubmodel,unitab)
      call hm_get_intv('c_FUNC', C_FUNC_ID, is_available,lsubmodel)
      call hm_get_floatv('cSCALE', Cscale, is_available,lsubmodel,unitab)

      call hm_get_floatv('RHO_TMD', rho_tmd, is_available,lsubmodel,unitab)
      call hm_get_floatv('C_SOLID', c_solid, is_available,lsubmodel,unitab)
      call hm_get_floatv('LAW5_PSH', psh, is_available,lsubmodel,unitab)

      call hm_get_intv('PLASEXP', PLASEXP, is_available,lsubmodel)
      call hm_get_intv('IFORM', IFORM, is_available,lsubmodel)

      !dimension / unit
      call hm_get_floatv_dim('RHO_TMD',density_unit,is_available,lsubmodel,unitab)

      rhor = pm(1)
      rhoi = pm(89)

      if(rho_tmd <= zero)then
         call ancmsg(MSGID=67,MSGTYPE=msgerror,ANMODE=aninfo,I1=uid,C1='/EOS/COMPACTION_TAB',C2='RHO_TMD MUST BE POSITIVE')
      endif

      if(c_solid <= zero)then
         call ancmsg(MSGID=67,MSGTYPE=msgerror,ANMODE=aninfo,I1=uid,C1='/EOS/COMPACTION_TAB',C2='C_SOLID MUST BE POSITIVE')
      endif

      if(PSCALE == zero) PSCALE = one

      if(cSCALE == zero) cSCALE = one

      !generate user function in data structure
      eos_param%ntable = 2
      allocate (eos_param%table(eos_param%ntable))
      eos_param%table(1)%notable = P_FUNC_ID
      eos_param%table(2)%notable = C_FUNC_ID
      x1scale   = one * density_unit
      x2scale   = one * density_unit
      x3scale   = one
      x4scale   = one
      x3vect(1) = zero
      x4vect(1) = zero
      fscale(1) = PSCALE
      fscale(2) = CSCALE
      call eos_table_copy(eos_param ,x2vect   ,x3vect   ,x4vect   ,         &
                x1scale  ,x2scale  ,x3scale  ,x4scale  ,fscale   ,         &
                ntable   ,table    ,ierror   ,uid)

      ! chekc pysical meaning : slope at the last point of the compaction curve vs max unloading modulus
      !NDIM = 1
      NPT = size(eos_param%table(1)%X(1)%VALUES)
      X1 = eos_param%table(1)%X(1)%VALUES(NPT-1)
      Y1 = eos_param%table(1)%Y1D(NPT-1)
      X2 = eos_param%table(1)%X(1)%VALUES(NPT)
      Y2 = eos_param%table(1)%Y1D(NPT)
      DX = X2-X1
      DY = Y2-Y1
      DYDX = ZERO
      IF (X1 /= X2) THEN
        DYDX = DY/DX
      ENDIF
      slope_end = DYDX
      if(slope_end > c_solid**2)then
        call ancmsg(MSGID=67,MSGTYPE=msgerror,ANMODE=aninfo,I1=uid,C1='/EOS/COMPACTION_TAB', &
                    C2='UNLOAD MODULUS C_SOLID**2 TOO LOW REGARDING COMPACTION CURVE')
      end if
      ! we must also check the slope at the intersection point with line of slope c_solid**2 crossing (rho_tmd,0)
      ! see rhomax_plastic calculation below

      !Default values
      ! Iform : formulation flag for unload behavior
      if(iform /= 1)then
        iform=1 !default
      endif

      !plasexp check and default
      if(plasexp /= 0 .and. plasexp /= 1)then
        plasexp = 0
      end if

      if(pm(79)==zero)pm(79)=three100

      ! P0=P(rho0)
      xvec1(1:1,1) = rhoi
      call table_mat_vinterp(eos_param%table(1),1,1,vartmp(1,1),xvec1,yy,slope)
      P0 = yy(1)
      ! SSP0 = Cunload(rho_0)
      call table_mat_vinterp(eos_param%table(2),1,1,vartmp(1,2),xvec1,yy,slope)
      SSP0 = yy(1)

      ! RHOMAX_PLASTIC
      ! intersection between
      !  1) the unloading curve (slope c_solid² crossing (rho_tmd,0)) :  eqn : yy1 = Puser(xx)
      !  2) and the compaction curve : yy(rho)                           eqn : yy2 = c_solid**2*xx - c_solid**2*rho_tmd
      ! newton iteration :  solving f(xx) = yy1-yy2 = 0
      NITER = 10
      ITER = 0
      TOL = em12
      rho_ = half*(X1+X2) !middle of the last segment
      residu = ep20
      DO WHILE (ITER <= NITER .AND. RESIDU > TOL)
        xvec1(1:1,1) = rho_
        call table_mat_vinterp(eos_param%table(1),1,1,vartmp(1,1),xvec1,yy,slope)
        Puser = yy(1)
        ff = Puser - c_solid**2 * (rho_-rho_tmd)
        df = slope(1) - c_solid**2
        rho_ = rho_ - ff/df
        residu = abs(ff/df)
      END DO
      rhomax_plastic = rho_ ! intersection point
      print *, "rhomax_plastic = ", rhomax_plastic


      !integer parameters
      eos_param%nuparam = 5
      eos_param%niparam = 2
      eos_param%nfunc = 0
      eos_param%ntable = 2
      call eos_param%construct() !allocations

      !real parameters
      eos_param%uparam(1) = rho_tmd
      eos_param%uparam(2) = c_solid
      eos_param%uparam(3) = PSH
      eos_param%uparam(4) = RHOMAX_PLASTIC

      !integer parameters
      eos_param%iparam(1) = iform
      eos_param%iparam(2) = plasexp

      !ssp0
      pm(27) = max(ssp0, pm(27))
      pm(23) = zero ! e0
      pm(31) = p0-psh
      pm(37) = em20
      pm(104)= p0-psh

      write(iout,1000)
      if(is_encrypted)then
        WRITE(IOUT,'(5X,A,//)')'CONFIDENTIAL DATA'
      else
        write(iout,1500) RHO_TMD,C_SOLID,PLASEXP
        write(iout,2000) P_FUNC_ID, PSCALE
        write(iout,2500) c_FUNC_ID, cSCALE
        write(iout,3000) rhomax_plastic
      endif

      return
! ----------------------------------------------------------------------------------------------------------------------

 1000 format(&
      5X,'  COMPACTION TABULATED EOS    ',/,&
      5X,'  ------------------------    ',/)
 1500 format(&
      5X,'RHO_TMD . . . . . . . . . . . . . . . . =',E12.4/&
      5X,'C_SOLID . . . . . . . . . . . . . . . . =',E12.4/&
      5X,'PLASTIC EXPANSION FLAG . . . . . . . .  =',I10)
 2000 format( &
      5X,'PLASTIC COMPACTION CURVE RHO_C(P_C) :',/ &
      5X,'  FUNCTION ID . . . . . . . . . . . . . =',I10/ &
      5X,'  SCALE FACTOR  . . . . . . . . . . . . =',E12.4)
 2500 format( &
      5X,'ELASTIC UNLOADING CURVE RHO_UNLOAD(C_UNLOAD) :',/ &
      5X,'  FUNCTION ID . . . . . . . . . . . . . =',I10/ &
      5X,'  SCALE FACTOR  . . . . . . . . . . . . =',E12.4)
 3000 format( &
      5X,'COMPUTED VALUES FROM COMPACTION CURVE :',/ &
      5X,'  RHOMAX_PLASTIC. . . . . . . . . . . . =',E12.4/)

      end subroutine hm_read_eos_compaction_tab
! ----------------------------------------------------------------------------------------------------------------------

      end module hm_read_eos_compaction_tab_mod
