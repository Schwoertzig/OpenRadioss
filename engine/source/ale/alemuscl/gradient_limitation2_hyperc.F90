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
      module gradient_limitation2_hyperc_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Compressive HYPER-C method
!! \details
      subroutine gradient_LIMITATION2_hyperc(trimat, ale_connect, numelq, nft, nel)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use ale_mod , only : ale
      use precision_mod , only : wp
      use constant_mod , only : zero, em12, ep20, one
      use ale_connectivity_mod , only : t_ale_connectivity
      use alemuscl_mod
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer, intent(in) :: trimat !< phase id
      type(t_ale_connectivity), intent(in) :: ale_connect ! ale connectivity buffer (elem/elem)
      integer,intent(in) :: numelq ! total number of quads
      integer,intent(in) :: nft !group shifting (elems)
      integer,intent(in) :: nel ! number of elem in the current gorup
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      integer :: i, ii, jj
      integer :: node_id
      real(kind=WP) :: reduc_factor(trimat), nodal_reduc_factor, yn, zn, valnode
      real(kind=WP) :: grad2,grad3,norm2,delta_up,delta_down,dalpha,dy,dz,r,phi_hyper,epsilon,proj
      integer :: itrimat,iadj,kk
      real(kind=WP) :: yk, zk
      real(kind=WP) :: beta, beta_max
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
      !!! Limiting process for the computed gradient -> maximum principle
      !!! and stability purposes
      epsilon = EM12
      beta = ale%inter_capturing%beta
      beta_max = beta

      do i = 1,nel
         ii = i + nft
         !element centroid
         yk = alemuscl_buffer%elcenter(ii,2) ; zk = alemuscl_buffer%elcenter(ii,3)
         reduc_factor = ep20
         do itrimat = 1, trimat
            if(abs(alemuscl_buffer%grad(ii,2,itrimat)) + abs(alemuscl_buffer%grad(ii,3,itrimat)) > zero) then
                   grad2 = alemuscl_buffer%grad(ii,2,itrimat)
                   grad3 = alemuscl_buffer%grad(ii,3,itrimat)
                   norm2 = grad2*grad2 + grad3*grad3
                   if (norm2 > epsilon) then
                    ! delta_up   : jump upwind <-> centroid
                    ! delta_down : jump centroid <-> downwind
                    ! using dotproduct (gradient direction)
                     delta_up   = zero
                     delta_down = zero
                     !--- loop over adjacent elem to identidy up/down
                     do kk = 1, 4
                       iadj = ale_connect%ee_connect%connected( ale_connect%ee_connect%iad_connect(ii) + kk - 1 )
                       if (iadj > 0 .and. iadj <= numelq) then
                         dalpha = alemuscl_buffer%volume_fraction(iadj,itrimat) - alemuscl_buffer%volume_fraction(ii,itrimat)
                         !  projecting along gradient direction
                         dy = alemuscl_buffer%elcenter(iadj,2) - yk
                         dz = alemuscl_buffer%elcenter(iadj,3) - zk
                         proj = dy*grad2 + dz*grad3
                         if (proj > zero) then
                           delta_down = max(delta_down, dalpha)
                         else if (proj < -zero) then
                           delta_up   = max(delta_up, -dalpha)
                         endif
                       endif
                     enddo
                     ! slope ratio r
                     r = zero
                     if (abs(delta_down) > epsilon) r = delta_up / (delta_down+em12)
                     ! limiteur hyper-c
                     ! phi_hyper = max(zero, min(beta*r, one))*beta_max
                     phi_hyper = min(beta, max(zero, beta*r))
                     ! applying limiter
                     reduc_factor(itrimat) = phi_hyper
                   else
                     reduc_factor(itrimat) = zero
                   end if
            else
               reduc_factor(itrimat) = zero
            endif
         enddo  ! itrimat = 1, trimat

         do itrimat = 1, trimat
            if(abs(alemuscl_buffer%grad(ii,2,itrimat)) + abs(alemuscl_buffer%grad(ii,3,itrimat)) > zero) then
               !!!   gradient limiter
               alemuscl_buffer%grad(ii,2,itrimat) = reduc_factor(itrimat) * alemuscl_buffer%grad(ii,2,itrimat)
               alemuscl_buffer%grad(ii,3,itrimat) = reduc_factor(itrimat) * alemuscl_buffer%grad(ii,3,itrimat)
            endif
         enddo
      enddo  ! i = 1,nel
! ----------------------------------------------------------------------------------------------------------------------
      end subroutine gradient_limitation2_hyperc
      end module gradient_limitation2_hyperc_mod
