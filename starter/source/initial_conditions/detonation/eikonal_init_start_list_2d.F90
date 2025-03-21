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
      module eikonal_init_start_list_2d_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief initialize narrow band for fast marching method depending on user input
!! \details
        subroutine eikonal_init_start_list_2d(nstart, start_elem_list, start_elem_tdet, detonators, numelq, numnod, &
                                              nod2elq, knod2elq, ale_connectivity, elem_list_bij, neldet, xel, x, &
                                              nummat, npropm, pm, nixq, ixq)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod, only : zero, ep20
          use detonators_mod , only : detonators_struct_
          use ale_connectivity_mod , only : t_ale_connectivity
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
          integer,intent(inout) :: nstart
          integer,intent(in) :: numelq, numnod !< array size
          integer,allocatable,dimension(:),intent(out) :: start_elem_list
          my_real,allocatable,dimension(:),intent(out) :: start_elem_tdet
          type (detonators_struct_),intent(in) :: detonators
          integer,intent(in) :: nod2elq(4*numelq) !< connectivity node->quad
          integer,intent(in) :: knod2elq(numnod+1) !< quad_id -> index in nod2elq (size 4)
          type (t_ale_connectivity), intent(inout) :: ale_connectivity
          integer,intent(in) :: elem_list_bij(numelq)
          integer,intent(in) :: neldet
          my_real,intent(in) :: xel(3,neldet)
          my_real,intent(in) :: x(3,numnod)
          integer,intent(in) :: npropm, nummat
          my_real,intent(in) :: pm(npropm,nummat)
          integer,intent(in) :: nixq
          integer,intent(in) :: ixq(nixq,numelq)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: idet, ndet_pts
          integer :: ii,jj,ielem
          integer :: start_id
          integer :: I_shadow_flag
          integer :: inod, nnod, nod_id
          integer :: num_adj, adjacent_elem(16) !max :4*4 (for a given node, attached quad : max4, Manathan neighborhood:max+8, diagonal neighborhood:max+4  (4+8+4=16)
          integer,allocatable,dimension(:) :: itag_elem
          integer iad1,lgth !< variable for ale_connectivity elem-elem buffer
          integer :: iev, iel, ie
          my_real :: dy,dz,dl,tdet,vel
          my_real :: ydet,zdet

          my_real,allocatable,dimension(:) :: tmp_tdet
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
        allocate(itag_elem(numelq)); !tag to check if elem was already found

        allocate(tmp_tdet(neldet))
        tmp_tdet(:) = ep20

        ! loop over detonation point with shadowing option (I_shadow_flag=1)
        ! and build ADJACENT NEIGHBORHOOD (3-STAGE-PROCESS)
        ndet_pts = detonators%n_det_point
        do idet = 1, ndet_pts
          I_shadow_flag = detonators%point(idet)%shadow
          nnod = detonators%point(idet)%shadow
          itag_elem(1:numelq) = 0
          ! loop over detonation points
          do inod=1,nnod
            nod_id = detonators%point(idet)%nodlist(inod)
            adjacent_elem(1:16) = 0
            num_adj = 0 !max 16
            ! --- FIRST STAGE
            ! loop over attached elems
            !ilen = knod2elq(nod_id+2)-knod2elq(nod_id+1)  !KNOD2ELQ(NQQ1)+1,KNOD2ELQ(NQQ1+1)
            do ii=KNOD2ELQ(nod_id)+1,KNOD2ELQ(nod_id+1)
              ielem = NOD2ELQ(ii)
              if(itag_elem(ielem) == 0)then
                itag_elem(ielem) = 100
                num_adj = num_adj + 1
                adjacent_elem(num_adj) = ielem
                print *, nod_id, ielem
              end if
            end do
            ! --- SECOND STAGE
            ! Manathan neighborhood
            do ii=1,num_adj
              ie = adjacent_elem(ii)
              iad1 = ale_connectivity%ee_connect%iad_connect(ie)
              lgth = ale_connectivity%ee_connect%iad_connect(ie+1) - iad1
              do jj=1,lgth
                iev = ale_connectivity%ee_connect%connected(iad1 + jj - 1)
                if(iev == 0)cycle
                iel = elem_list_bij(iev)
                if(itag_elem(iel) == 0)then
                   itag_elem(iel) = 10
                   num_adj = num_adj + 1
                   adjacent_elem(num_adj) = iel
                 end if
              end do
            end do
            ! --- THIRD STAGE
            ! Adjacent neighborhood
            do ii=1,num_adj
              ie = adjacent_elem(ii)
              if(itag_elem(ie) /= 10)cycle
              iad1 = ale_connectivity%ee_connect%iad_connect(ie)
              lgth = ale_connectivity%ee_connect%iad_connect(ie+1) - iad1
              do jj=1,lgth
                iev = ale_connectivity%ee_connect%connected(iad1 + jj - 1)
                if(iev == 0)cycle
                iel = elem_list_bij(iev)
                if(itag_elem(iel) < 10)then
                   itag_elem(iel) = itag_elem(iel) + 1
                 end if
              end do
            end do
            do ii=1,numelq
              if(itag_elem(ii) == 2)then
                 num_adj = num_adj + 1
                 adjacent_elem(num_adj) = ii
              endif
            end do

            ! ---INIT WITH RADIAL DISTANCE
            vel = zero
            do ii=1,num_adj
              iel = adjacent_elem(ii)
              vel = max(vel,pm(38,ixq(1,iel))) ! chapman jouget velocity
            end do

            ydet = x(2,nod_id)
            zdet = x(3,nod_id)
            do ii=1,num_adj
              iel = adjacent_elem(ii)
              iel = elem_list_bij(iel)
              dy = ydet - xel(1,iel)  !xel-storage is 1:y 2:z
              dz = zdet - xel(2,iel)
              dl = sqrt(dy*dy + dz*dz)
              tdet = dl /vel
              tmp_tdet(iel) = min (tmp_tdet(iel), tdet)
            end do

            ! --RESULT
            do ii=1,num_adj
              print *, "adjacent_elem(ii)=",adjacent_elem(ii), tmp_tdet(elem_list_bij(adjacent_elem(ii)))
            end do

            print *, " "
            ! recuperer les elements attaches a ce noeud
            ! ajouter les voisins (si mat_det valide)
            ! ajouter les voisins en diagonale (si mat_det valide)
            ! initialiser avec propagation radiale  (tdetdoit avoir ete initialisé à 1e20)
          enddo


        end do

        nstart = 0
        do ii=1,neldet
          if(tmp_tdet(ii) /= ep20)then
            nstart = nstart + 1
          end if
        end do

        allocate(start_elem_list(nstart))
        allocate(start_elem_tdet(nstart))

        nstart = 0
        do ii=1,neldet
          if(tmp_tdet(ii) /= ep20)then
            nstart = nstart + 1
            start_elem_list(nstart) = ii
            start_elem_tdet(nstart) = tmp_tdet(ii)
          end if
        end do

        deallocate(itag_elem)
        deallocate(tmp_tdet)

      end subroutine eikonal_init_start_list_2d
! ----------------------------------------------------------------------------------------------------------------------
      end module eikonal_init_start_list_2d_mod
