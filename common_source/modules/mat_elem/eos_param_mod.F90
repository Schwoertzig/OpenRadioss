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

      module eos_param_mod

! ======================================================================================================================
!! \brief module to define data structure for viscosity model parameters in materials
!! \details 

      use table4d_mod
      use names_and_titles_mod

      implicit none
!
#include "my_real.inc"
!
!=======================================================================      
      
      type eos_param_
        character(len=nchartitle) :: title = ''  !< eos model input name
        integer     :: nuparam                   !< number of real value paraameters
        integer     :: niparam                   !< number of int value parameters
        integer     :: nuvar                     !< number of internal state variables
        integer     :: nfunc                     !< number of local functions in material
        integer     :: ntable                    !< number of local function tables
        
        my_real        ,dimension(:) ,allocatable :: uparam  !< real value eos parameter table
        integer        ,dimension(:) ,allocatable :: iparam  !< int  value eos parameter table
        integer        ,dimension(:) ,allocatable :: func    !< function table in eos models
        type(table_4d_),dimension(:) ,allocatable :: table   !< local function tables
      
      end type eos_param_   
!
!---------------
      end module eos_param_mod
