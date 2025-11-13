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
      module alemain_cutcell_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Main Subroutine called by the resol loop to treat ALE MULTIMAT with cut-cell interface tracking
!! \details ...
        subroutine  alemain_cutcell(nixq,numelq,ixq, numnod, x, ale_connect, ncycle, &
                                    ityptstt, neltstt, t1s, tt, &
                                    ngroup, ngrnod, igrnod, nparg, iparg, dt_scale, dt1, dt2t, multi_cutcell)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod, only : pi,zero,one,em01               !module containing all constant insitialized with either single or double precision
          use precision_mod, only : wp                            !provides kind for eigther single or double precision (wp means working precision)
          use ale_connectivity_mod , only : t_ale_connectivity    !data structure for ale elem-elem connectivities
          use groupdef_mod , only : group_
          use ale_mod , only : ALE
          use debug_mod , only : ITAB_DEBUG
          use multi_cutcell_mod, only : multi_cutcell_struct, allocate_multi_cutcell_type
          use multi_cutcell_solver_mod, only : initialize_solver_multicutcell, update_fluid_multicutcell
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
          ! input
          integer,intent(in) :: nixq                                    !< array size for ixq array
          integer,intent(in) :: numelq                                  !< number of quad in the input file
          integer,intent(in) :: ixq(nixq,numelq)                        !< elem connectivity (1:7 : mat_id,n1,n2,n3,n4,pid,user_id for each elem 1:numelq)
          integer,intent(in) :: numnod                                  !< number of nodes in the input file
          integer,intent(in) :: ngroup                                  !< number of groups (elem buffer data structure)
          integer,intent(in) :: nparg                                   !< array size for IPARG data structure
          integer,intent(in) :: iparg(nparg,ngroup)                     !< parameter for groups (used with elbuf)
          integer,intent(in) :: ncycle                                  !< resol cycle number (time loop)
          real(kind=wp),intent(in) :: x(3,numnod)                       !< node coordinates
          real(kind=wp),intent(in) :: tt                                !< current time
          type(t_ale_connectivity), intent(inout) :: ale_connect
          integer,intent(in) :: ngrnod                                  !< number of group of nodes(array size for igrnod)
          type(group_)  ,dimension(ngrnod)  :: igrnod                   !< group of nodes (data structure)
          real(kind=WP), intent(in) :: dt_scale
          real(kind=wp),intent(in) :: dt1

          !output
          real(kind=wp),intent(out) :: dt2t
          type(multi_cutcell_struct), intent(inout) :: multi_cutcell        !< element buffer (storage for output files : pressure, density, velocity, ...)
          integer,intent(inout) :: ityptstt                             !< element type imposing the time step (here 2)
          integer,intent(inout) :: neltstt                              !< element user id imposing the time step
          real(kind=WP),intent(inout) :: t1s                            !< simulation time
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: part_id, polyg_id !< user parameters
          integer :: nel      !< number of element in the group
          integer :: nft      !< shift (group partitionning nel<=128)
          integer :: ng       !< loop on groups
          integer :: ii       !< local id in the group
          integer :: elem_iid !< internal id
          integer :: elem_uid !< user id
          integer :: iad2,lgth, iadj, JJ !< ale_connectivity usage (elem-elem)
          real(kind=WP) :: result_density !< temp
          real(kind=WP) :: gamma
          integer :: sign !Not used for now
          integer :: nb_phase !TODO should be input
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Precondition
! ----------------------------------------------------------------------------------------------------------------------
          if(ALE%solver%multimat%is_defined_mmale3 == 0) return  ! otherwise interface tracking with cut-cell is required
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
! [ the code must be indented with 2 spaces]
! [ the code must be commented ]
! [ routines must be short, and should not exceed 200 lines for leaf routines, 1000 lines for main routines]

! [ separators can be used between blocks of code]

                 ! ---------------------------DO PRE-TRETMENT HERE
                 ! -----------------------------------------------
                 !user input
                 part_id = ALE%solver%multimat%list(1)%part_id ! by default we can take into account all elem and ignore part_id
                 polyg_id = ALE%solver%multimat%list(1)%surf_id

                 ! EXAMPLE OF USAGE
                 if(ncycle == 0)then
                   ! elem connecvitivy
                   print *, "elem iid=1 : uid,n1,n2,n3,n4 =,", 1, ixq(7, 1), ixq(2:5, 1)   !iid : internal id, uid : user id
                   print *, "elem iid=1  uid(n1),uid(n2),uid(n3),uid(n4) =,", ITAB_DEBUG(ixq(2:5, 1))   !iid : internal id, uid : user id

                   ! polygon
                   print *, "user polygon internal identifier : ", polyg_id
                   print *, "     number of nodes : ", igrnod(polyg_id)%nentity
                   print *, "     list of nodes   : ", igrnod(polyg_id)%entity (1:igrnod(polyg_id)%nentity)

                   !ale_connectivity (elem-elem)
                   ELEM_IID=1
                   IAD2 = ALE_CONNECT%ee_connect%iad_connect(ELEM_IID)
                   LGTH = ALE_CONNECT%ee_connect%iad_connect(ELEM_IID+1) - IAD2
                   DO JJ=1,LGTH
                     IADJ = ALE_CONNECT%ee_connect%connected(IAD2 + JJ - 1)
                     IF(IADJ > 0)THEN
                         !there is an adjacent elem on face J, its elem_iid is IADJ
                      END IF
                   ENDDO

                   ! node coordinates
                   print *, "node iid=3, x,y,z=", X(1:3, 3), "uid=", ITAB_DEBUG(3)

                 end if ! ncycle == 0


                 ! -------------------------------CALL 2D POC HERE
                 ! -----------------------------------------------
                 gamma = 1.4
                 sign = 1 !not used for now

                 if (dt1>0) then
                  call update_fluid_multicutcell(n2d, numelq, numeltg, numnod, ixq, ixtg, x, itab, ALE_CONNECT, &
                         multi_cutcell%grid, multi_cutcell%phase_vely, multi_cutcell%phase_velz, multi_cutcell%phase_rho, multi_cutcell%phase_pres, &
                         gamma, dt1, dt_scale, sign, dt2t, multi_cutcell%sound_speed)
                 else
                  !Initialization
                  nb_phase = 2
                  call allocate_multi_cutcell_type(nb_phase, numelq + numeltg, multi_cutcell)
                  call initialize_solver_multicutcell(n2d, numelq, numeltg, numnod, ixq, ixtg, x, itab, ALE_CONNECT, &
                              ALE%solver%multimat%nb, ALE%solver%multimat%list(:)%surf_id, igrnod, multi_cutcell%grid)
                 end if
                   




                 ! -------------------------DO POST-TREATMENT HERE
                 ! -----------------------------------------------
                 ! required for time loop
                 dt2t = em01  ! time step used computed or imposed by the poc
                 neltstt = 1  ! element user id imposing the minimal time step
                 ityptstt = 2 !element type imposing the minimal time step
                 T1S = TT     ! required to go to next cycle

                 ! example of USAGE
                 result_density = em01
                 do ng=1,ngroup
                   nel = iparg(2,ng)
                   nft = iparg(3,ng) ! shift
                   do ii=1,nel
                     elem_iid = ii+nft
                     elbuf(ng)%gbuf%rho(ii) = result_density ! affect density for elem_iid here
                   end do
                 enddo



! ----------------------------------------------------------------------------------------------------------------------
        end subroutine alemain_cutcell
      end module alemain_cutcell_mod
