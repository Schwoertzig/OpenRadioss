#ifndef UPDATE_CLIPPED_H
#define UPDATE_CLIPPED_H

#include "my_real.h"
#include "Polygon2D.h"
#include "Polyhedron3D.h"
#include "array_int.h"
#include <stdint.h>

Polyhedron3D* build_space2D_time_cell(const Polygon2D *fn, const my_real *vsx, const my_real* vsy, uint64_t size_vs, const my_real dt, bool split, Array_int *list_del_pts);
void update_solid(Polygon2D **solid, Polyhedron3D** solid3D, const my_real* vec_move_solidx, const my_real* vec_move_solidy, my_real dt, \
                    my_real minimal_length, my_real maximal_length, my_real minimal_angle);
#endif