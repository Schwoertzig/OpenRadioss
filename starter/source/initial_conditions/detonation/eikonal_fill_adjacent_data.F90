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
      module eikonal_fill_adjacent_data_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
        subroutine eikonal_fill_adjacent_data(ie, ALE_CONNECTIVITY,neldet, &
                                              tdet, tdet_adj, vel,vel_adj, xel,xel_adj, numel,elem_list_bij, &
                                              updown,num_new_activated, list_new_activated,  updown_adj)
!! \brief In order to apply Godunov operator the adjacent data must be used. This subroutine is collecting these data (positions, velocities, arrival times)
!! \details ALE_EE_CONNECT is used to retrive data from adjacent element (ALE,EULER only)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use ale_connectivity_mod , only : t_ale_connectivity
          use constant_mod , only : ep20
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
          integer,intent(in) :: ie
          type (t_ale_connectivity), intent(inout) :: ale_connectivity
          integer,intent(in) :: neldet
          my_real,intent(inout) :: tdet_adj(6)
          my_real,intent(inout) :: xel_adj(3,6)
          my_real,intent(inout) :: vel_adj(6)
          integer,intent(inout) :: updown_adj(6)
          integer,intent(in) :: numel
          integer, intent(in) :: elem_list_bij(numel)  ! size
          integer,intent(inout) :: updown(neldet)
          integer,intent(inout) :: num_new_activated, list_new_activated(6)
          my_real,intent(in) :: vel(neldet),xel(3,neldet)
          my_real,intent(inout) :: tdet(neldet)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: jj
          integer :: iad1,lgth
          integer iel
          integer iev
          my_real :: s,dx,dy,dl
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
           num_new_activated = 0
           list_new_activated(1:6)=0
           iad1 = ale_connectivity%ee_connect%iad_connect(ie)
           lgth = ale_connectivity%ee_connect%iad_connect(ie+1) - iad1
           do jj=1,lgth
             tdet_adj(jj) = ep20
             vel_adj(jj) = vel(ie)
             updown_adj(jj) = -1
             iev = ale_connectivity%ee_connect%connected(iad1 + jj - 1)
             if(iev == 0)cycle
             iel = elem_list_bij(iev)
             tdet_adj(jj) = tdet(iel)
             vel_adj(jj) = vel(iel)
             xel_adj(1:2,jj) = xel(1:2,iel)
             updown_adj(jj) = updown(iel)
             if(updown(iel) == -1) then
               num_new_activated = num_new_activated + 1
               list_new_activated(num_new_activated) = iel
               updown(iel) = 0
               s=max(vel(iel),vel(ie))
               s=1./s
               dx = xel(1,ie)-xel_adj(1,jj)
               dy = xel(2,ie)-xel_adj(2,jj)
               dl = sqrt(dx*dx+dy*dy)
               tdet(iel) = tdet(ie) + dl*s
             end if
           ENDDO
        end subroutine eikonal_fill_adjacent_data
! ----------------------------------------------------------------------------------------------------------------------
        end module eikonal_fill_adjacent_data_mod
