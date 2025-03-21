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
      module eikonal_solver_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Solver for the Eikonal equation (arrival time depending on medium velocity) using the FAST MARCHING METHOD
        subroutine eikonal_solver(ixq      , nixq     , numelq  , &
                                  ixtg     , nixtg    , numeltg , &
                                  ixs      , nixs     , numels  , &
                                  x        , numnod   , ipri    , &
                                  elbuf_tab, ngroup   , nparg   , &
                                  nod2eltg , knod2eltg, &
                                  nod2elq  , knod2elq , &
                                  nod2els  , knod2els , &
                                  iparg    , ale_connectivity, npropm, nummat, pm, n2d, detonators)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use eikonal_solver_2d_mod
          use elbufdef_mod, only : elbuf_struct_
          use ale_connectivity_mod , only : t_ale_connectivity
          use constant_mod , only : zero, ep20
          use detonators_mod , only : detonators_struct_
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
#include "units_c.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer,intent(in) :: ngroup  !< number of groups
          integer,intent(in) :: numnod  !< number of nodes
          integer,intent(in) :: nixq    !< array size ixq
          integer,intent(in) :: numelq  !< number of quad elems (solid)
          integer,intent(in) :: nixs    !< array size ixs
          integer,intent(in) :: numeltg !< number of tria elems (solid)
          integer,intent(in) :: nixtg   !< array size ixtg
          integer,intent(in) :: numels  !< number of hexa elems (solid)
          integer,intent(in),target :: ixq(nixq,numelq) !< quad connectivities
          integer,intent(in),target :: ixtg(nixtg,numeltg) !< quad connectivities
          integer,intent(in),target :: ixs(nixs,numels) !< quad connectivities
          integer,intent(in) :: nparg !< array size
          integer,intent(in) :: iparg(nparg,ngroup)
          my_real,intent(in) :: x(3,numnod) !< node coordinates
          integer,intent(in) :: npropm, nummat
          my_real,intent(in) :: pm(npropm,nummat)
          type (elbuf_struct_), target, dimension(ngroup) :: elbuf_tab
          type (t_ale_connectivity), intent(inout) :: ale_connectivity
          integer,intent(in) :: n2d
          integer,intent(in) :: ipri
          type (detonators_struct_),intent(in) :: detonators
          integer,intent(in) :: nod2eltg(3*numeltg)
          integer,intent(in) :: nod2elq(3*numelq)
          integer,intent(in) :: nod2els(3*numels)
          integer,intent(in) :: knod2eltg(numnod+1)
          integer,intent(in) :: knod2elq(numnod+1)
          integer,intent(in) :: knod2els(numnod+1)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: i,nel,ng,nft,mpr,iel
          my_real :: tdet
          integer,pointer,dimension(:,:) :: ix
          integer :: nix
          integer :: I_shadow_flag
          integer :: idet, mat_det
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Precondition
! ----------------------------------------------------------------------------------------------------------------------
          if( detonators%n_det == 0 )return              !NO DETONATORS => RETURN
          if( numels + numelq + numeltg*n2d == 0)return   !NO SOLID ELEMS => RETURN
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          if(nod2els(1)==knod2els(1))then
            ! unused argument used
          end if
          if(nod2eltg(1)==knod2eltg(1))then
            ! unused argument used
          end if


          do idet = 1, detonators%n_det_point

            I_shadow_flag = detonators%point(idet)%shadow
            If(I_shadow_flag == 0)cycle
            Mat_det = detonators%point(idet)%mat

            ! FAST MARCHING METHOD
            if(n2d == 0) then
              !call eikonal_solver_3d(ixs,nixs,numels,x, numnod,&
              !                       elbuf_tab,ngroup,nparg,iparg,ale_connectivity,npropm,nummat,pm,&
              !                       detonators, idet, nod2els, knod2els,&
              !                       numels, numnod)
            else
              call eikonal_solver_2d(ixq,nixq,numelq,x,numnod, &
                                     elbuf_tab,ngroup,nparg,iparg,ale_connectivity,npropm,nummat,pm,&
                                     detonators, idet, nod2elq, knod2elq)
            end if

          end do

       !---------------------------------!
       !    UNINT ELEMS (NO DETONATOR)   !
       !---------------------------------!
       do ng = 1,ngroup
         nel = iparg(2,ng)
         if(elbuf_tab(ng)%gbuf%g_tb > 0)then
           do i=1,nel
             tdet = elbuf_tab(ng)%gbuf%tb(i)
             if(tdet == -ep20)then
               elbuf_tab(ng)%gbuf%tb(i) = zero
             end if
           enddo
         endif
       end do

       !---------------------------------!
       !            PRINTOUT             !
       !---------------------------------!
       if(n2d == 0)then
          ix => ixs(1:nixs,1:numels)
          nix = nixs
        elseif(numelq > 0)then
          ix => ixq(1:nixq,1:numelq)
          nix = nixq
        else
          ix => ixtg(1:nixtg,1:numeltg)
          nix = nixq
        end if

       if(ipri >= 3)then
         mpr =0
        do ng = 1,ngroup
          nel = iparg(2,ng)
          nft = iparg(3,ng)
          do i=1,nel
            mpr = mpr+1
            iel = ix(nix,i+nft)
            if(mpr == 1) write(iout,500)
            tdet=-elbuf_tab(ng)%gbuf%tb(i)
            write(iout,510) nel,tdet
            if(mpr == 50) mpr=0
          end do
        end do
       endif

! ----------------------------------------------------------------------------------------------------------------------
500  FORMAT(//, &
       5X, 'DETONATION TIMES FOR JWL ELEMENTS' /, &
       5X, '---------------------------------' //, &
       5X, 'ELEMENT DETONATION TIME' /)
  510 FORMAT(5X,I10,E15.5)
! ----------------------------------------------------------------------------------------------------------------------

      end subroutine eikonal_solver
! ----------------------------------------------------------------------------------------------------------------------

      end module eikonal_solver_mod
