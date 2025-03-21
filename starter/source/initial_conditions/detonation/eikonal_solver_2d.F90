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
      module eikonal_solver_2d_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Solver for the Eikonal equation (arrival time depending on medium velocity)
!! \details updown(:) array is used to define narrow band 1:frozen, 0:narrow_band, -1:far
        subroutine eikonal_solver_2d(ixq,nixq,numelq,x,numnod, &
                                     elbuf_tab,ngroup,nparg,iparg,ale_connectivity,npropm,nummat,pm,&
                                     detonators, idet, nod2elq, knod2elq)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod, only : zero, ep20, fourth
          use elbufdef_mod, only : elbuf_struct_
          use ale_connectivity_mod , only : t_ale_connectivity
          use eikonal_sort_narrow_band_mod , only : eikonal_sort_narrow_band
          use eikonal_remove_first_mod , only : eikonal_remove_first
          use eikonal_fill_adjacent_data_mod , only : eikonal_fill_adjacent_data
          use eikonal_init_start_list_2d_mod , only : eikonal_init_start_list_2d
          use detonators_mod , only : detonators_struct_
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
          integer,intent(in) :: ngroup !< number of groups
          integer,intent(in) :: nixq   !< array size ixq
          integer,intent(in) :: numnod !< number of nodes
          integer,intent(in) :: numelq !< number of quad elems (solid)
          integer,intent(in) :: ixq(nixq,numelq) !< quad connectivities
          integer,intent(in) :: nparg !< array size
          integer,intent(in) :: iparg(nparg,ngroup)
          my_real,intent(in) :: x(3,numnod) !< node coordinates
          integer,intent(in) :: npropm, nummat
          my_real,intent(in) :: pm(npropm,nummat)
          type (elbuf_struct_), target, dimension(ngroup) :: elbuf_tab
          type (t_ale_connectivity), intent(inout) :: ale_connectivity
          type (detonators_struct_),intent(in) :: detonators
          integer,intent(in) :: idet
          integer,intent(in) :: nod2elq(3*numelq)
          integer,intent(in) :: knod2elq(numnod+1)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: ng
          integer :: ity
          integer :: mlw
          integer :: nel
          integer :: neldet
          integer :: nft
          integer :: mat_det
          integer :: nod_id
          integer :: start_id
          integer :: iel
          integer :: i  !<local elem id
          integer :: ii,jj !< loops
          integer :: iad1,LGTH,ie, iev !< local variable for adjacent cell loop
          integer :: inod(8) !tria 3<=8,  quad 4<=8, hexa 8
          integer :: n_queue
          integer, allocatable, dimension(:) :: idx_ng  !< group identifier for a given solid elem
          integer, allocatable, dimension(:) :: idx_i   !< local id in group
          integer, allocatable, dimension(:) :: elem_list
          my_real, allocatable, dimension(:) :: vel
          my_real :: vel_adj(6)   ! tria 3<6, quad:4<6, hexa 6
          my_real, allocatable, dimension(:) :: tdet   !< detonation time
          my_real :: tdet_adj(6)  ! tria 3<6, quad:4<6, hexa 6
          my_real, allocatable, dimension(:,:) :: Xel   !< centroids coordinates
          my_real :: xel_adj(3,6) ! tria 3<6, quad:4<6, hexa 6
          integer :: updown_adj(6)
          integer, allocatable, dimension(:) :: priority_queue_id
          my_real, allocatable, dimension(:) :: priority_queue_tt
          integer,allocatable, dimension(:) :: updown ! -1 down, 1:up, 0:narrow_band
          integer :: nstart !< number of deotnation points (centroids)  ! can be adapt later from mesh nodes to elem centroid (spherical wave from node to centroid)
          integer,allocatable,dimension(:) :: start_elem_list
          my_real,allocatable,dimension(:) :: start_elem_tdet
          integer, allocatable, dimension(:) :: elem_list_bij
          integer num_new_activated, list_new_activated(6)
          my_real :: dx,dy,dl,s, tmp
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

         mat_det = detonators%point(idet)%mat

         !numbering
         neldet = 0
         do ng=1,ngroup
           mlw = iparg(1,ng)
           nel = iparg(2,ng)
           nft = iparg(3,ng)
           ity = iparg(5,ng)
           if(ity == 2)then
             if(mlw == 5 .or. mlw == 97)then
               !if(ixq(1,nft+1) == mat_det .or. mat_det == 0)then
                 neldet = neldet + nel
              end if
           endif
         enddo

         ! size neldet
         allocate(elem_list(neldet)) ; elem_list(:) = 0
         allocate(idx_ng(neldet))    ; idx_ng(:) = 0
         allocate(idx_i(neldet))     ; idx_i(:) = 0
         allocate(updown(neldet))    ; updown(:) = -1
         allocate(tdet(neldet))      ; tdet(:) = ep20
         allocate(vel(neldet))       ; vel(:) = ep20 !lowest chosen for stability
         allocate(Xel(3,neldet))     ; Xel(:,:) = zero
         allocate(priority_queue_id(neldet)) ; priority_queue_id(:) = 0
         allocate(priority_queue_tt(neldet)) ; priority_queue_tt(:) = ep20

         ! size numelq
         allocate(elem_list_bij(numelq)) ; elem_list_bij = 0

         ! list of relevant centroids
         !   and group id and local id in this group : elem_list(k) -> (idx_ng(k), idx_i(k))   ! used to set burning time ELBUF_TAB(NG)MGBUF%TB(I)
         neldet = 0
         do ng=1,ngroup
           mlw = iparg(1,ng)
           nel = iparg(2,ng)
           nft = iparg(3,ng)
           ity = iparg(5,ng)
           if(ity == 2)then
             if (mlw == 5 .or. mlw == 97)then
               !if(ixq(1,nft+1) /= mat_det .and. mat_det /= 0)cycle
               do i=1,nel
                 !building list
                 neldet = neldet + 1
                 elem_list(neldet) = i+nft       ! elem_list     : [1..neldet] -> [1..numelq]
                 idx_ng(neldet) = ng
                 idx_i(neldet) = i
                 elem_list_bij(i+nft) = neldet   ! elem_list_inv : [1..numelq] -> [1..neldet]
                 !centroid coordinates
                 inod(1:4) = ixq(2:5, i+nft)
                 xel(1, neldet) = fourth * sum(x(2,inod(1:4))) ! x-center
                 xel(2, neldet) = fourth * sum(x(3,inod(1:4))) ! y-center
                 xel(3, neldet) = zero
                 !medium velocity
                 vel(neldet) = pm(38,ixq(1,i+nft)) ! chapman jouget velocity
               enddo
              endif
           endif
         enddo

         call eikonal_init_start_list_2d(nstart, start_elem_list, start_elem_tdet, detonators, numelq, numnod, &
                                         nod2elq, knod2elq, ale_connectivity, elem_list_bij, neldet, xel, x,&
                                         nummat,npropm,pm, nixq, ixq)


        print *,"  ci dessous skip traitement si mat_det pas coherent"

         ! initial detonation locations
         ! mark first upwind points (centroids)
         do jj=1,nstart
           updown(start_elem_list(jj)) = 1
           tdet(start_elem_list(jj)) = start_elem_tdet(jj)
         enddo ! next jj

         ! initial narrow band
         !    mark first points in the narrow band (close)
         !    by the way list their detonation time
         n_queue = 0
         do ii=1,neldet
           if(updown(ii) == 1)then
             !tag adjacent cells with updown -1
             ie = elem_list(ii)
             IAD1 = ALE_CONNECTIVITY%ee_connect%iad_connect(IE)
             LGTH = ALE_CONNECTIVITY%ee_connect%iad_connect(IE+1) - IAD1
             DO JJ=1,LGTH
               IEV = ALE_CONNECTIVITY%ee_connect%connected(IAD1 + JJ - 1)
               IF(IEV == 0)CYCLE
               iel = elem_list_bij(iev)
               if(updown(iel) == -1) then
                 updown(iel) = 0
                 n_queue = n_queue + 1
                 s=max(vel(iel),vel(ie))
                 s=1./s
                 dx = xel(1,ie)-xel(1,iel)
                 dy = xel(2,ie)-xel(2,iel)
                 dl = sqrt(dx*dx+dy*dy)
                 priority_queue_id(n_queue) = iel
                 priority_queue_tt(n_queue) = tdet(ie) + dl*s
                 tdet(iel) = priority_queue_tt(n_queue)
               end if
             ENDDO
           end if
         enddo ! next ii
         call eikonal_sort_narrow_band(priority_queue_id,priority_queue_tt,n_queue)

         ! main loop -------------------------------------------------------------------------
         do while (n_queue > 0)

           ie = priority_queue_id(1) ! list of  priority_queue_tt is already sorted

           !fill adjacent data (arrival time and velocity)
           ! and awake adjacent which were 'far'
           call eikonal_fill_adjacent_data(ie, ALE_CONNECTIVITY,neldet, &
                                           tdet,tdet_adj,vel,vel_adj,xel,xel_adj,numelq,elem_list_bij, &
                                           updown, num_new_activated, list_new_activated,  updown_adj)

           !update arrival time on this point(centroid)
           call eikonal_Godunov_Operator_2d_bis(xel(1,ie), tdet(ie), xel_adj, tdet_adj, 4, Vel(ie), Vel_adj, &
                                            updown(ie))

           priority_queue_tt(1) = tdet(ie)
           call eikonal_remove_first(priority_queue_id,priority_queue_tt,n_queue)

           do ii=1,num_new_activated
             n_queue = n_queue + 1
             ie = list_new_activated(ii)
             priority_queue_id(n_queue) = ie
             priority_queue_tt(n_queue) = tdet(ie)
           end do

         enddo !wend
         ! end of main loop -------------------------------------------------------------------------

         ! initialize element buffer (arrival times)
         do ii=1,neldet
           ng = idx_ng(ii)
           i  = idx_i(ii)
           tmp = -elbuf_tab(ng)%gbuf%tb(i)
           tmp = min (tmp, tdet(ii))
           elbuf_tab(ng)%gbuf%tb(i) = -tmp
         end do

         !DEALLOCATE
         deallocate(start_elem_list)
         deallocate(start_elem_tdet)
         deallocate(elem_list)
         deallocate(idx_ng)
         deallocate(idx_i)
         deallocate(updown)
         deallocate(tdet)
         deallocate(vel)
         deallocate(Xel)
         deallocate(priority_queue_id)
         deallocate(priority_queue_tt)
         deallocate(elem_list_bij)

        end subroutine eikonal_solver_2d
! ----------------------------------------------------------------------------------------------------------------------



! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Godunov operator for fast marching method in 2D using finite differences
!! \details Operator valid only for structured meshes (alignes adjacent points)
        subroutine eikonal_Godunov_Operator_2d(xel, tt, xel_adj, tt_adj, n_adj, Velocity, Velocity_adj, updown)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod , only : zero, one, ep20, two, four, em06
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
          integer, intent(in) :: n_adj ! number of adjacent points
          integer,intent(inout) :: updown
          my_real, intent(in) :: Xel(3) ! centroids coordinates
          my_real, intent(in) :: Xel_adj(3, n_adj) ! centroids coordinates of adjacent points
          my_real, intent(inout) :: tt ! arrival time
          my_real, intent(in) :: tt_adj(n_adj) ! arrival time of adjacent points
          my_real, intent(in) :: Velocity ! velocity on current point
          my_real, intent(in) :: Velocity_adj(n_adj) ! velocity on adjacent points

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: k,l
          my_real :: dx, dy
          my_real :: tt_candidate
          my_real :: a, b
          my_real :: s
          my_real :: aa,bb,cc,delta
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          dx=xel(1)-xel_adj(1,1)
          dx=0.1
          dy=0.1
          ! Initialize gradients
          a = min(tt_adj(1),tt_adj(3))
          b = min(tt_adj(2),tt_adj(4))
          if (a == ep20 .and. b == ep20) return
          if (a == ep20)then
            s = max(Velocity_adj(2) , Velocity_adj(4))
            s = max(Velocity, s)
            s = one / s
            tt_candidate = b+s*dx
            tt = min(tt, tt_candidate)
          elseif (b == ep20) then
            s = max(Velocity_adj(1) , Velocity_adj(3))
            s = max(Velocity, s)
            s = one / s
            tt_candidate = a+s*dx
            tt = min(tt, tt_candidate)
          else
            s = maxval(Velocity_adj)
            s= max(s, Velocity)
            s = one / s
            k=1
            if(a==tt_adj(3))k=3  !supposing x axis fior adj elems 1,3
            l=2
            if(b==tt_adj(4))l=4  !supposing x axis fior adj elems 2,4
            aa = dx*dx+dy*dy
            bb = -two*(dy*dy*a+dx*dx*b)
            cc = dy*dy*a*a + dx*dx*b*b - dx*dx*dy*dy*s*s
            delta = bb*bb - 4.*aa*cc
            if(delta >=0)then
              tt_candidate = (-bb + sqrt(delta)) / two / AA
            else
              tt_candidate = (-bb + zero ) / two / AA
            end if
            tt = min(tt, tt_candidate)
          end if

          updown = 1

        end subroutine eikonal_Godunov_Operator_2d
! ----------------------------------------------------------------------------------------------------------------------



! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Godunov operator for fast marching method in 2D using gradient reconstruction (unstructured mesh)
!! \details Finite differences cannot be used since adjacent points are not necessarily aligned.
        subroutine eikonal_Godunov_Operator_2d_bis(xel, tt, xel_adj, tt_adj, n_adj, Velocity, Velocity_adj, updown)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod , only : zero, one, ep20, two, four, em06
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
          integer, intent(in) :: n_adj ! number of adjacent points
          integer,intent(inout) :: updown
          my_real, intent(in) :: Xel(3) ! centroids coordinates
          my_real, intent(in) :: Xel_adj(3, n_adj) ! centroids coordinates of adjacent points
          my_real, intent(inout) :: tt ! arrival time
          my_real, intent(in) :: tt_adj(n_adj) ! arrival time of adjacent points
          my_real, intent(in) :: Velocity ! velocity on current point
          my_real, intent(in) :: Velocity_adj(n_adj) ! velocity on adjacent points

! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: k,l
          my_real :: dx, dy, dist
          my_real :: tt_candidate
          my_real :: a, b
          my_real :: s
          my_real :: delta
          my_real :: A1,A2,B1,B2,C1,C2,AA,BB,CC,DENOM
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          ! Initialize gradients
          a = min(tt_adj(1),tt_adj(3))
          b = min(tt_adj(2),tt_adj(4))
          if (a == ep20 .and. b == ep20) return
          if (a == ep20)then
            l=2
            if(b==tt_adj(4))l=4  !supposing x axis fior adj elems 2,4
            dx = (xel(1)-xel_adj(1,l))
            dy = (xel(2)-xel_adj(2,l))
            dist = sqrt(dx*dx+dy*dy)
            s = max(Velocity_adj(2) , Velocity_adj(4))
            s = max(Velocity, s)
            s = one / s
            tt_candidate = b+s*dist
            tt = min(tt, tt_candidate)
          elseif (b == ep20) then
            k=1
            if(a==tt_adj(3))k=3  !supposing x axis fior adj elems 1,3
            dx = (xel(1)-xel_adj(1,k))
            dy = (xel(2)-xel_adj(2,k))
            dist = sqrt(dx*dx+dy*dy)
            s = max(Velocity_adj(1) , Velocity_adj(3))
            s = max(Velocity, s)
            s = one / s
            tt_candidate = a+s*dist
            tt = min(tt, tt_candidate)
          else
            s = maxval(Velocity_adj)
            s= max(s, Velocity)
            s = one / s
            k=1
            if(a==tt_adj(3))k=3  !supposing x axis fior adj elems 1,3
            l=2
            if(b==tt_adj(4))l=4  !supposing x axis fior adj elems 2,4

            A1 = (xel_adj(2,k)-xel_adj(2,l)) ;
            A2 = (xel_adj(1,l)-xel_adj(1,k)) ;
            B1 = two*( (xel_adj(2,l)-xel(2))*tt_adj(k) - (xel_adj(2,k)-xel(2))*tt_adj(l) )*A1
            B2 = two*( (xel_adj(1,k)-xel(1))*tt_adj(l) - (xel_adj(1,l)-xel(1))*tt_adj(k) )*A2
            C1 = ( (xel_adj(2,l)-xel(2))*tt_adj(k) - (xel_adj(2,k)-xel(2))*tt_adj(l) )
            C2 = ( (xel_adj(1,k)-xel(1))*tt_adj(l) - (xel_adj(1,l)-xel(1))*tt_adj(k) )

            DENOM =  (xel_adj(1,k)-xel(1))*(xel_adj(2,l)-xel(2)) - (xel_adj(1,l)-xel(1))*(xel_adj(2,k)-xel(2))
            DENOM = DENOM*DENOM

            A1 = A1*A1
            A2 = A2*A2
            C1 = C1*C1
            C2 = C2*C2

            AA = (A1+A2)/DENOM
            BB = (B1+B2)/DENOM
            CC = (C1+C2)/DENOM - s*s

            delta = BB*BB-FOUR*AA*CC

            if(delta >= zero)then
              tt_candidate = (-BB + sqrt(delta)) / two / AA
            else
              tt_candidate = (-BB + ZERO) / TWO / AA
            end if
            tt = min(tt, tt_candidate)

          end if

          updown = 1

        end subroutine eikonal_Godunov_Operator_2d_bis
! ----------------------------------------------------------------------------------------------------------------------
      end module eikonal_solver_2d_mod
