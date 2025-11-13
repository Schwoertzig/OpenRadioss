#include "my_real.inc"

module polygon_mod
type Point2D
  SEQUENCE 
  my_real :: y
  my_real :: z
end type Point2D

type Point3D
  SEQUENCE 
  my_real :: y
  my_real :: z
  my_real :: t
end type Point3D
end module polygon_mod
