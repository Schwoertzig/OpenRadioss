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
      subroutine ebcs11(nseg,iseg,segvar, &
                        a,v,w,x, &
                        liste,nod,irect,ielem,iface, &
                        la,ms,stifn,ebcs,iparg,elbuf_tab,ixq,ixs,ixtg, &
                        fsavsurf,lsky,fsky,elem_adress,time,iparit,dt1, &
                        numels, numelq, numeltg, numnod, nparg, ngroup, nixs, nixq, nixtg, nsurf, iale, n2d)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   modules
! ----------------------------------------------------------------------------------------------------------------------
      use ebcs_mod
      use elbufdef_mod
      use multi_fvm_mod
      use segvar_mod
      use ale_mod , only : ale
      use th_surf_mod , only : th_surf_num_channel
      use constant_mod , only : zero, one, third, fourth, half, two
! ----------------------------------------------------------------------------------------------------------------------
!                                                   implicit none
! ----------------------------------------------------------------------------------------------------------------------
      implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   include files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer,intent(in) :: n2d !< 2d/3d flag
      integer,intent(in) :: iale !< ale flag
      integer,intent(in) :: numeltg, numels, numelq, numnod, nparg, ngroup, nixs, nixq, nixtg, nsurf !< array sizes
      my_real, intent(in) :: dt1 !< time step
      my_real, intent(in) :: time !< simulation time
      integer,intent(in) :: iparit !< /parith/on flag
      integer, dimension(4,nseg), intent(in) :: elem_adress ! adress for fsky array (only used with parith/on)
      integer,intent(in) :: lsky !< array size
      my_real, dimension(8,lsky), intent(inout) :: fsky ! acceleration array for parith/on option
      my_real, intent(inout) :: fsavsurf(th_surf_num_channel,nsurf)
      integer,intent(in) :: nseg,nod,iseg(nseg),liste(nod),irect(4,nseg),ielem(nseg),iface(nseg)
      integer,intent(in) :: ixq(nixq,numelq),ixs(nixs,numels),ixtg(nixtg,numeltg)
      my_real,intent(inout) :: ms(numnod)
      my_real,intent(inout) :: a(3,numnod),v(3,numnod),w(3,numnod),x(3,numnod),la(3,nod),stifn(numnod)
      type(t_ebcs_propergol), intent(inout) :: ebcs
      integer :: iparg(nparg,ngroup)
      type(elbuf_struct_), target, dimension(ngroup) :: elbuf_tab
      type(t_segvar),intent(inout) :: segvar
! ----------------------------------------------------------------------------------------------------------------------
!                                                   local variables
! ----------------------------------------------------------------------------------------------------------------------
      type(g_bufel_), pointer :: gbuf
      type(l_bufel_)  ,pointer :: lbuf
      integer ii,is,kk,k1,ll,kseg,nn(4),nng(4),num,kty,klt,mft,ngrp,iloc,npt,ivoi,idx(6),ix(4)
      integer icf_2d(2,4), icf_3d(4,6), jj, isubmat, ipos, mtn
      integer isolnod
      my_real orient,rho,c,fc,roc,alp,vol,mass,mass_face,mass_face_old
      my_real x13,y13,z13,x24,y24,z24,nx,ny,z,s
      my_real roou,enou,vmx,vmy,vmz,fluxi,fluxo,p,dvx,dvy,dvz
      my_real v0(3,nod)
      my_real xsn,ysn,zsn
      my_real xn, yn, zn, vold, vnew, pold, padj
      my_real dp0,mach,pp,ssp,rhoc2,vel_front,fac1,fac2
      my_real param_a, param_n, param_q, surf,phase_alpha(21),phase_rho(21), phase_eint(21)
      my_real face_force
      my_real :: param_rho0s  !< initial propergol density
      my_real :: dmass_g     !< mass increment (burnt propergol)
      my_real :: dvol_s      !< burnt volume
      my_real :: deint_g     !< released energy from burnt propergol
      my_real :: tmp(3)

      logical bfound
      type(buf_mat_)  ,pointer :: mbuf
      integer :: adress  !< adress for parit/on array
      integer :: param_surf_id !< ebcs surface identifier
      integer :: nbsubmat

      data icf_2d  /1,2,2,3,3,4,4,1/
      data icf_3d  /1,4,3,2,3,4,8,7,5,6,7,8,1,2,6,5,2,3,7,6,1,5,8,4/

      INTEGER N0PHAS
      PARAMETER (N0PHAS = 04)

      INTEGER NVPHAS
      PARAMETER (NVPHAS = 23)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   body
! ----------------------------------------------------------------------------------------------------------------------

      !--- user parameters
      param_a = ebcs%a
      param_n = ebcs%n
      param_q = ebcs%q
      param_rho0s = ebcs%rho0s
      param_surf_id = ebcs%surf_id

      !--- nodal velocities at face nodes
      do ii=1,nod
        num=liste(ii)
        v0(1,ii)=v(1,num)
        v0(2,ii)=v(2,num)
        v0(3,ii)=v(3,num)
      enddo
      if(iale == 1)then
        do ii=1,nod
          num=liste(ii)
          v0(1,ii)=v0(1,ii)-w(1,num)
          v0(2,ii)=v0(1,ii)-w(2,num)
          v0(3,ii)=v0(1,ii)-w(3,num)
        enddo
      endif

      !--- reset working array for internal forces
      do ii=1,nod
        num=liste(ii)
        la(1,ii)=zero
        la(2,ii)=zero
        la(3,ii)=zero
      enddo


      do is=1,nseg
        kseg=abs(iseg(is))
        orient=float(iseg(is)/kseg)
        !---outward normal
        if(n2d == 0)then
          jj = iface(is)
          ix(1)=ixs(icf_3d(1,jj)+1,ielem(is))
          ix(2)=ixs(icf_3d(2,jj)+1,ielem(is))
          ix(3)=ixs(icf_3d(3,jj)+1,ielem(is))
          ix(4)=ixs(icf_3d(4,jj)+1,ielem(is))
          x13=x(1,ix(3))-x(1,ix(1))
          y13=x(2,ix(3))-x(2,ix(1))
          z13=x(3,ix(3))-x(3,ix(1))
          x24=x(1,ix(4))-x(1,ix(2))
          y24=x(2,ix(4))-x(2,ix(2))
          z24=x(3,ix(4))-x(3,ix(2))
          xn=y13*z24-z13*y24
          yn=z13*x24-x13*z24
          zn=x13*y24-y13*x24
          fac2=one/sqrt(xn*xn+yn*yn+zn*zn)
          xn = xn*fac2
          yn = yn*fac2
          zn = zn*fac2
          surf = half/fac2
          if(ix(4) == ix(3))then ; npt=3;fac1=third; else; npt=4;fac1=fourth; endif
        else
          fac1=half
          npt=2
          jj = iface(is)
          if(numeltg>0)then
            ix(1)  = ixtg(icf_2d(1,jj)+1,ielem(is))
            ix(2)  = ixtg(icf_2d(2,jj)+1,ielem(is))
          else
            ix(1)  = ixq(icf_2d(1,jj)+1,ielem(is))
            ix(2)  = ixq(icf_2d(2,jj)+1,ielem(is))
          endif
          xn = zero
          yn = -(-x(3,ix(2))+x(3,ix(1)))
          zn =  (-x(2,ix(2))+x(2,ix(1)))
          fac2 = one/sqrt(yn*yn+zn*zn)
          yn=yn*fac2
          zn=zn*fac2
          surf = one/fac2
        endif
        tmp(1:3) = zero
        nn(1)=irect(1,is)
        nn(2)=irect(2,is)
        nn(3)=irect(3,is)
        nn(4)=irect(4,is)
        nng(1:4)=liste(nn(1:4))
        do kk=1,npt
          tmp(1) = tmp(1) + v0(1,nng(kk))
          tmp(2) = tmp(2) + v0(2,nng(kk))
          tmp(3) = tmp(3) + v0(3,nng(kk))
        enddo
        vnew = fac1 * (tmp(1)*xn + tmp(2)*yn + tmp(3)*zn)
        !-- adjacent state
        bfound = .false.
        ivoi = ielem(is)
        do ngrp=1,ngroup
          kty = iparg(5,ngrp)
          klt = iparg(2,ngrp)
          mft = iparg(3,ngrp)
          isolnod = iparg(28,ngrp)
          if(n2d==0)then
            if(kty /= 1)cycle
          else
            if(kty == 2)then
              isolnod=4
            elseif(kty == 7)then
              isolnod=3
            else
              cycle
            endif
          endif
           if (ivoi <= klt+mft)then
             bfound = .true.
             exit
           endif
        enddo
        if(.not.bfound)cycle !next i
        gbuf => elbuf_tab(ngrp)%gbuf
        lbuf => elbuf_tab(ngrp)%bufly(1)%lbuf(1,1,1)
        mtn = iparg(1,ngrp)
        !adjacent pressure
        iloc = ivoi-mft-1
        do kk=1,6
          idx(kk) = klt*(kk-1)
        enddo
        padj = -third*(gbuf%sig(idx(1)+iloc+1) + gbuf%sig(idx(2)+iloc+1) + gbuf%sig(idx(3)+iloc+1))
        !density
        rho = gbuf%rho(iloc+1)
        vol = gbuf%vol(iloc+1)
        mass = vol * rho
        !adjacent sound speed
        ssp = lbuf%ssp(iloc+1)
        !volume fracions and submat state
        if(mtn == 51)then
          mbuf => elbuf_tab(ngrp)%bufly(1)%mat(1,1,1)
          do isubmat=1,4
            ipos = 1
            kk = (n0phas + (isubmat-1)*nvphas +ipos-1) * klt  +  iloc+1
            phase_alpha(isubmat) = mbuf%var(kk)
            ipos = 9
            kk = (n0phas + (isubmat-1)*nvphas +ipos-1) * klt  +  iloc+1
            phase_rho(isubmat) = mbuf%var(kk)
            ipos = 8
            kk = (n0phas + (isubmat-1)*nvphas +ipos-1) * klt  +  iloc+1
            phase_eint(isubmat) = mbuf%var(kk)
          enddo!next itrimat
          vol = vol * phase_alpha(1) ! burnt gas is supposed to be submat 1
          mass = vol * phase_rho(1)
        endif
        !-- formulation
        !burnt propergol volume
        pp = padj
         !adjacent elem
        if(pp <= zero)then
          vel_front = zero
          dvol_s = zero
          dmass_g = zero
        else
          vel_front = param_a*exp(param_n*log(pp))
          dvol_s = surf*vel_front*dt1
          dmass_g = dvol_s * param_rho0s
        endif
        vold = ebcs%vold(is)
        ebcs%vold(is) = vel_front
        if(time == zero) then
            vold = zero
        endif
        ! ghost cell update ( upwind/aconve() )
        segvar%rho(kseg)=param_rho0s
        segvar%eint(kseg)=param_q*param_rho0s
        ! burnt gas is supposed to be submat 1
        if(segvar%nbmat > 1)then
          nbsubmat = segvar%nbmat
          phase_rho(1) = param_rho0s
          segvar%phase_rho(2:nbsubmat,kseg)=phase_rho(2:nbsubmat)
          phase_eint(1) = param_q*param_rho0s
          segvar%phase_eint(2:nbsubmat,kseg)=phase_eint(2:nbsubmat)
          segvar%phase_alpha(1:nbsubmat,kseg)=phase_alpha(1:nbsubmat)
        endif
        !-- expand pressure to face nodes
        face_force = pp*surf                                        !pp for equilibrium
        !mass = mass + dmass_g
        mass_face = mass*(npt*one)/(isolnod*one)
        !vold = zero
        if(dt1 > zero)face_force = face_force + (vel_front) * mass_face / dt1   !pp= pp+ dp to input corresponding propergol
        !expand pressure loading to segment nodes
        !face_force = pp*surf
        do kk=1,npt
          la(1,nn(kk)) = la(1,nn(kk)) - fac1* (face_force) * xn
          la(2,nn(kk)) = la(2,nn(kk)) - fac1* (face_force) * yn
          la(3,nn(kk)) = la(3,nn(kk)) - fac1* (face_force) * zn
          if (ms(nn(kk)) > zero)then
            stifn(nn(kk))=stifn(nn(kk))+(two*(surf*rho*ssp)**2)/ms(nn(kk))
          endif
        enddo
         !/th/surf (massflow, velocity, pressure)
         fsavsurf(2,param_surf_id) = fsavsurf(2,param_surf_id) + rho*surf*vel_front     !rho.s.u = dm/dt
         fsavsurf(3,param_surf_id) = fsavsurf(3,param_surf_id) + surf*vel_front         !s.u
         fsavsurf(4,param_surf_id) = fsavsurf(4,param_surf_id) + surf*padj              !s.p
         fsavsurf(6,param_surf_id) = fsavsurf(6,param_surf_id) + rho*surf*vel_front*dt1 ! m<-m+dm (cumulative)
         ! -------------
         ! for parith/on option : update forces in fsky array (specific assembly)
         if(iparit > 0) then
           do kk=1,npt
             adress = elem_adress(kk,is) ! adress of fsky array for element is and node kk
             fsky(1,adress) = -face_force*xn*fac1
             fsky(2,adress) = -face_force*yn*fac1
             fsky(3,adress) = -face_force*zn*fac1
             fsky(4:8,adress) = zero
           enddo
         endif
         ! -------------
      enddo

      ! numerical staggered scheme : vitesse.F ;  v[n+1] = v[n] + int(acc, t=t[n],t[n+1] )   ;   acc = F/m
      !   to impose v[n+1] = Vel_Front
      !           -->  v[n] is reset in velocity subroutine for nodes related to /EBCS/PROPERGOL

      ! -------------
      ! for parith/off option : update directly the acceleration array a() : no specific assembly
      if(iparit == 0) then
        do ii=1,nod
           num=liste(ii)
           a(1,num)=a(1,num)+la(1,ii)
           a(2,num)=a(2,num)+la(2,ii)
           a(3,num)=a(3,num)+la(3,ii)
        enddo
      endif
      ! -------------


      return
      end
