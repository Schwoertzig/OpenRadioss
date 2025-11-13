#include "my_real.inc"

module grid2D_struct_multicutcell_mod
  use polygon_mod

  type grid2D_struct_multicutcell
    logical                     :: close_cells !Determine whether a cell is close to the interface
    logical                     :: is_narrowband
    !my_real, dimension(4)       :: length_edge
    my_real, dimension(4)       :: lambda_per_edge
    my_real                     :: lambdan_per_cell
    my_real                     :: lambdan_per_cell_target
    my_real                     :: lambdanp1_per_cell
    my_real                     :: lambdanp1_per_cell_target
    type(Point2D)               :: normal_intern_face_space
    my_real                     :: normal_intern_face_time
    type(Point2D)               :: p_normal_intern_face_space
    my_real                     :: p_normal_intern_face_time
    my_real                     :: area
  end type grid2D_struct_multicutcell

  type ConservativeState2D
    my_real :: rho
    my_real :: rhoE
    my_real :: rhovz
    my_real :: rhovy 
  end type ConservativeState2D

  type ConservativeFlux2D
    my_real, dimension(4) :: rho
    my_real, dimension(4) :: rhoE
    my_real, dimension(4) :: rhovz
    my_real, dimension(4) :: rhovy 
  end type ConservativeFlux2D

  contains 
  function compute_area(N,X) result(area)  !Area enclosed by polygon P, by the "shoelace" method.
      implicit none
      integer :: N !The number of points.
      my_real :: X(2, N) !The points.
      my_real :: area

      area = sum(X(1,1:N - 1)*X(2,2:N)) + X(1,N)*X(2,1) &
             - sum(X(1,2:N)*X(2,1:N - 1)) - X(1,1)*X(2,N)
      area = abs(area/2.0)  !The midpoint formula requires a halving.
  end function compute_area

  function makegrid(NUMELQ, NUMELTG, NUMNOD, IXQ, IXTG, X, i_cell) result(grid)
    implicit none
    integer :: N2D, NUMELQ, NUMELTG, NUMNOD, i_cell
    integer, dimension(:,:) :: IXQ, IXTG
    my_real, dimension(:,:) :: X
    type(grid2D_struct_multicutcell) :: grid

    grid%close_cells = .true.
    grid%is_narrowband = .true.

    if (NUMELQ > 0) then
      grid%area = compute_area(4, X(2:3, IXQ(2:5, i_cell)))
    else
      grid%area = compute_area(3, X(2:3, IXTG(2:4, i_cell)))
    end if

  end function makegrid

end module grid2D_struct_multicutcell_mod
