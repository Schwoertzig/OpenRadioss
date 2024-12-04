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
!! \brief Reset velocity on boundary surface when /EBCS/PROPERGOL is used.
!! \details The general formulation v[n+1] = v[n] + acc.dt.  To impose v[n+1] = vel_burning_front, v[n] is set to 0
!! \details   it leads to v[n+1] = 0 + acc.dt = vel_burning_front. acc is computed from internal force which is computed in ebcs subroutine
        subroutine ebcs_reset_vel(ebcs_tab, v, numnod)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use ebcs_mod, only : t_ebcs_tab, t_ebcs_propergol, t_ebcs
          use constant_mod , only : zero
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
          type(t_ebcs_tab), target, intent(inout) :: ebcs_tab
          my_real, intent(inout) :: v(3,numnod)
          integer,intent(in) :: numnod
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: i,ii !< loops
          integer :: typ  !< ebcs type : 11 is /EBCS/PROPERGOL
          integer :: node_id  !< node internal identifier
          class(t_ebcs), pointer :: ebcs
          class(t_ebcs_propergol), pointer :: twf
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          do i = 1, ebcs_tab%nebcs
            ebcs => ebcs_tab%tab(i)%poly
            typ = ebcs%type
            if (typ == 11) then
              select type (ebcs)
                type is (t_ebcs_propergol)
                  twf => ebcs
                  !--- nodal velocities at face nodes
                  do ii = 1, twf%nb_node
                    node_id = twf%node_list(ii)
                    v(1,node_id) = zero
                    v(2,node_id) = zero
                    v(3,node_id) = zero
                  enddo
              end select
            endif
          enddo
! ----------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------
        end subroutine ebcs_reset_vel
