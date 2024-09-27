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
      subroutine ebcs_main_post(segvar,a,v,w,x,ms,stifn,iparg,&
                           elbuf_tab,ebcs_tab,ixq,ixs,ixtg,&
                           fsky,fsavsurf,time,dt1,numnod,nparg,ngroup, &
                           nixq, nixs, nixtg, numels, numelq, numeltg,&
                           lsky,nsurf,iparit,iale,n2d)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
      use elbufdef_mod
      use groupdef_mod
      use ebcs_mod
      use th_surf_mod , only : th_surf_num_channel
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit None
! ----------------------------------------------------------------------------------------------------------------------
       implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
      integer,intent(in) :: n2d !< 2d/3d fla
      integer,intent(in) :: iale !< ale flag
      integer,intent(in) :: iparit !< parith/on flag
      integer,intent(in) :: lsky, nsurf !< array sizes
      integer,intent(in) :: nixq, nixs, nixtg, numels, numelq, numeltg  !< array sizes
      integer,intent(in) :: nparg, ngroup
      integer,intent(in) :: numnod
      my_real,intent(in) :: dt1 !time step
      my_real,intent(in) :: time !simulation time
      my_real,intent(inout) :: fsavsurf(th_surf_num_channel,nsurf)
      integer iparg(nparg,ngroup)
      integer,intent(in) :: ixq(nixq,numelq),ixs(nixs,numels),ixtg(nixtg,numeltg)
      my_real segvar(*),v(3,numnod),w(3,numnod),a(3,numnod),x(3,numnod),ms(numnod),stifn(numnod)
      type(elbuf_struct_), dimension(ngroup) :: elbuf_tab
      type(t_ebcs_tab), target, intent(inout) :: ebcs_tab
      my_real, dimension(8,lsky), intent(inout) :: fsky ! acceleration array for parith/on option
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local Variables
! ----------------------------------------------------------------------------------------------------------------------
      integer i,typ,isu,nseg,nod,j
      class(t_ebcs), pointer :: ebcs
      logical has_th
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

      do i = 1, ebcs_tab%nebcs
         if(.not.ebcs_tab%need_to_compute(i)) cycle
         ebcs => ebcs_tab%tab(i)%poly
         if(ebcs%is_multifluid)return
         typ = ebcs%type
         isu = ebcs%surf_id
         nseg = ebcs%nb_elem
         nod = ebcs%nb_node
         if (typ == 11) then
            select type (twf => ebcs)
             type is (t_ebcs_propergol)
              call ebcs11(nseg,twf%iseg,segvar, &
                          a,v,w,x,   &
                          twf%node_list,nod,twf%elem_list,twf%ielem,twf%iface, &
                          twf%la,ms,stifn,twf,iparg,elbuf_tab,ixq,ixs,ixtg,  &
                          fsavsurf,lsky, fsky, ebcs_parithon(i)%elem_adress,time,iparit,dt1, &
                          numels, numelq, numeltg,numnod, nparg, ngroup, nixs, nixq, nixtg, nsurf, iale, n2d)
            end select
         endif
      enddo

      return
      end subroutine ebcs_main_post
      
