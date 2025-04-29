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
      !||    eos_compaction_tab_mod   ../common_source/eos/eos_compaction_tab.F90
      !||--- called by ------------------------------------------------------
      !||    eosmain           ../common_source/eos/eosmain.F
      !||====================================================================
      module compaction_tab_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief  This subroutine contains numerical solving of eos_compaction_tab EOS
!----------------------------------------------------------------------------
!! \details STAGGERED SCHEME IS EXECUTED IN TWO PASSES IN EOSMAIN : IFLG=0 THEN IFLG=1
!! \details COLLOCATED SCHEME IS DOING A SINGLE PASS : IFLG=2
!! \details
!! \details  STAGGERED SCHEME
!! \details     EOSMAIN / IFLG = 0 : DERIVATIVE CALCULATION FOR SOUND SPEED ESTIMATION c[n+1] REQUIRED FOR PSEUDO-VISCOSITY (DPDE:partial derivative, DPDM:total derivative)
!! \details     MQVISCB            : PSEUDO-VISCOSITY Q[n+1]
!! \details     MEINT              : INTERNAL ENERGY INTEGRATION FOR E[n+1] : FIRST PART USING P[n], Q[n], and Q[n+1] CONTRIBUTIONS
!! \details     EOSMAIN / IFLG = 1 : UPDATE P[n+1], T[N+1]
!! \details                          INTERNAL ENERGY INTEGRATION FOR E[n+1] : LAST PART USING P[n+1] CONTRIBUTION
!! \details                            (second order integration dE = -P.dV where P = 0.5(P[n+1] + P[n]) )
!! \details  COLLOCATED SCHEME
!! \details     EOSMAIN / IFLG = 2 : SINGLE PASS FOR P[n+1] AND DERIVATIVES
!----------------------------------------------------------------------------
      !||====================================================================
      !||    eos_compaction_tab     ../common_source/eos/eos_compaction_tab.F90
      !||--- called by ------------------------------------------------------
      !||    eosmain         ../common_source/eos/eosmain.F
      !||--- uses       -----------------------------------------------------
      !||    constant_mod    ../common_source/modules/constant_mod.F
      !||    eos_param_mod   ../common_source/modules/mat_elem/eos_param_mod.F90
      !||====================================================================
      subroutine compaction_tab(&
                            iflag , nel   , pm     , off    , eint   , &
                            dvol  , mat   , psh    , dt1    , rho    , rho0  , &
                            pnew  , dpdm  , dpde   , rho_bak, &
                            npropm, nummat, nvareos, vareos , nvartmp, vartmp, &
                            eos_param)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
       use constant_mod , only : zero, em10, half, one, two, ep20
       use eos_param_mod , only : eos_param_
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
      integer,intent(in) :: nvareos
      my_real,intent(inout) :: vareos(nel,nvareos)
      integer,intent(in) :: nel !< number of element in the currenbt group
      integer,intent(in) :: npropm, nummat !< array sizes
      integer,intent(in) :: mat(nel), iflag
      my_real,intent(inout) :: pm(npropm,nummat) !< material data (real parameters)
      my_real,intent(inout) :: off(nel),eint(nel),dvol(nel)
      my_real,intent(inout) :: pnew(nel),dpdm(nel),dpde(nel)
      type(eos_param_),intent(in) :: eos_param !< data structure for EoS parameters
      my_real,intent(inout) :: rho_bak(nel) !< backup of mu for unloading
      my_real,intent(in) :: dt1 !< time step
      my_real,intent(in) :: rho(nel)  !< current density
      my_real,intent(in) :: rho0(nel) !< initial density
      integer ,intent(in) :: nvartmp                       !< size for vartmp
      integer ,dimension(nel,nvartmp) ,intent(inout) :: vartmp    !< vartmp is the history of index position on the user curve (optimization in order not to loop from first point at each cycle)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      integer i,mx,iter,niter
      my_real :: residu, ff, df
      my_real :: p0,psh(nel)
      my_real :: p(nel)

      my_real :: rho_tmd
      my_real :: c_solid

      my_real :: lambda

      !integer parameters
      integer :: iform
      integer :: plasexp

      my_real :: tol
      my_real :: pmin
      my_real :: rhomax_plastic

      my_real, dimension(nel,1) :: xvec1 !<temporary array for table interpolation

      my_real :: cunl(nel), c_prime(nel)
      my_real :: Puser(nel)
      my_real :: dpdr(nel)

      my_real :: b(nel)

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
       mx             = mat(1)
       psh(1:nel)     = pm(88,mx)
       p0             = pm(31,mx)
       pmin           = pm(37,mx)
       rho_tmd        = eos_param%uparam(1)
       c_solid        = eos_param%uparam(2)
       psh            = eos_param%uparam(3)
       rhomax_plastic = eos_param%uparam(4)

      !integer parameters
      iform = eos_param%iparam(1)
      plasexp = eos_param%iparam(2)

      IF(DT1 == zero)then
         call compaction_tab_init(eos_param,nel,nvartmp,vartmp,nvareos,vareos,rho,rho_tmd)
      ENDIF!(DT1 == zero)then

       if(iflag == 0)then

         do i=1,nel
           ! --- path along compaction curve
           if(rho(i) > rho_bak(i) .and. rho(i) < rhomax_plastic)then

               ! p(i) = P(rho)
               xvec1(1:1,1) = rho(i)
               call table_mat_vinterp(eos_param%table(1),1,1,vartmp(1,1),xvec1,Puser,dPdr)
               P(i) = Puser(1)
               P(i) = max(P(i), Pmin)

               ! update lambda
               if(iform == 1)then  !linear unloading
                 ! Updated lambda   (find lambda such as unloading path interesct compaction curve
                 iter = 0
                 residu = ep20
                 lambda=vareos(i,1)
                 do while(iter <= niter .and. residu > tol)
                   xvec1(1:1,1) = lambda
                   call table_mat_vinterp(eos_param%table(2),1,1,vartmp(1,2),xvec1,cunl,c_prime) ! get C_unload(labmda)
                   FF = P(i) - cunl(1)*cunl(1) * (rho(i)-lambda) - Pmin
                   DF = -two * cunl(1) * c_prime(1) * (rho(i)-lambda) + cunl(1)*cunl(1)
                   LAMBDA = LAMBDA - FF / DF
                   residu = abs(FF)/RHO_TMD
                   iter = iter + 1
                 enddo
               endif

               vareos(I,1) = min(lambda, rho_tmd)       ! lambda
               vareos(I,2) = cunl(1)*cunl(1)            ! c(lambda)**2
               dpdm(i) = rho0(i)*dPdr(i)                ! total derivative for sound speed
               rho_bak(i) = min(rhomax_plastic,rho(i))  ! history rho_plastic
           else

             ! path along unloading line
             if(iform == 1)then  !linear unloading
               if(rho(i) > rhomax_plastic)then
                 vareos(i,1) = rho_tmd
                 vareos(i,2) = c_solid*c_solid
               end if
               lambda = vareos(1,i)
               b(i) = vareos(i,2)
               P(i) = b(i)*(rho(i)-lambda) + Pmin
               P(i) = max(P(i), Pmin)
               dpdm(i) = rho0(i)*b(i)
             endif
           endif
         enddo !next i
         dpde(1:nel) = zero
         pnew(1:nel) = p(1:nel)*off(1:nel)-psh(1:nel)
         VAREOS(1:nel,3) = Pnew(1:nel)
      elseif(iflag == 1) THEN
        !----------------------------------------------------------------!
        !  FRACTURE  - MU_BAK                                            !
        !----------------------------------------------------------------!
        do i=1,nel
          pnew(i) = vareos(i,3) ! pressure already calculated (no energy increment since dpde=0)
          eint(i) = eint(i) - half*dvol(i)*(pnew(i)+psh(i) )
        enddo !next i
      endif

!------------------------
      return
      end subroutine compaction_tab
!------------------------




      subroutine compaction_tab_init(eos_param,nel,nvartmp,vartmp,nvareos,vareos,rho,rho_tmd)
!! \brief  This subroutine initialize VAREOS aray for initial time 0.0
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
       use constant_mod , only : zero, em10, half, one, two, ep20
       use eos_param_mod , only : eos_param_
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
      integer,intent(in) :: nvareos
      my_real,intent(inout) :: vareos(nel,nvareos)
      integer,intent(in) :: nel !< number of element in the currenbt group
      integer ,intent(in) :: nvartmp                       !< size for vartmp
      integer ,dimension(nel,nvartmp) ,intent(inout) :: vartmp    !< vartmp is the history of index position on the user curve (optimization in order not to loop from first point at each cycle)
      type(eos_param_),intent(in) :: eos_param !< data structure for EoS parameters
      my_real,intent(in) :: rho(nel),rho_tmd
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variable
! ----------------------------------------------------------------------------------------------------------------------
      my_real :: Puser(nel),dPdr(nel),Cunl(nel),slope(nel)
      my_real :: residu,lambda,ff,df,tol
      my_real, dimension(nel,1) :: xvec1 !<temporary array for table interpolation
      integer :: i, iter,niter
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
         niter = 10
         tol = em10

         !The current unloading speed-of-sound velocity
         xvec1(1:nel,1) = rho(1:nel)
         call table_mat_vinterp(eos_param%table(2),nel,nel,vartmp(1,2),xvec1,cunl,slope)
         VAREOS(1:nel,2) = Cunl(1:nel)*Cunl(1:nel)

         !solve the zero-crossing of the current unloading
         call table_mat_vinterp(eos_param%table(1),nel,nel,vartmp(1,1),xvec1,Puser,dPdr) !The current compaction pressure
          DO I=1,NEL
            IF(Puser(i) == zero)THEN
              VAREOS(I,1) = RHO(I) ! LAMBDA
            ELSE
              iter = 0
              residu = ep20
              lambda = rho(i)
              do while(iter <= niter .and. residu > tol)
                xvec1(1:1,1) = lambda
                call table_mat_vinterp(eos_param%table(1),1,1,vartmp(1,1),xvec1,Puser,dPdr) ! provides cunl=c(x) and c_prime=c'(x)
                FF = Puser(1)
                DF = dPdr(1)
                LAMBDA = LAMBDA - FF / DF
                residu = abs(FF)/RHO_TMD
                iter = iter + 1
              enddo
              VAREOS(I,1) = lambda
            END IF
          ENDDO

          VAREOS(1:nel,3) = Puser(1:nel)

      end subroutine compaction_tab_init


      end module compaction_tab_mod
