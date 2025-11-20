#ifndef BRIDGE_C_FORTRAN_H
#define BRIDGE_C_FORTRAN_H

#include "my_real.h"
#include "Polygon2D.h"


void launch_grb_();
void end_grb_();

void build_grid_from_points_fortran_(const my_real* x_v, const my_real* y_v, const long long int *signed_nb_pts);
void build_clipped_from_pts_fortran_(const my_real* x_v, const my_real* y_v, const long long int* limits, const long long *nb_polygons);

/// @brief 
/// @param grid 
/// @param clipped 
/// @param vec_move_solid 
/// @param size_vec_move_solid 
/// @param dt 
/// @param ptr_lambdas_arr 
/// @param height_lambdas_arr 
/// @param width_lambdas_arr 
/// @param ptr_big_lambda_n 
/// @param size_big_lambda_n 
/// @param ptr_big_lambda_np1 
/// @param size_big_lambda_np1 
/// @param mean_normal 
/// @param is_narrowband 
void compute_lambdas2d_fortran_(const my_real *dt, \
                        my_real *ptr_lambdas_arr,    \
                        my_real *ptr_big_lambda_n,   \
                        my_real *ptr_big_lambda_np1, \
                        Point3D *mean_normal, bool *is_narrowband);
void nb_pts_clipped_fortran_(long long int* signed_nb_pts_solid);
void nb_edge_clipped_fortran_(long long int* signed_nb_edges_solid);
void compute_normals_clipped_fortran_(my_real* normalVecx, my_real* normalVecy, \
                                        my_real* normalVecEdgex, my_real* normalVecEdgey, \
                                        my_real* min_pos_Se);
void smooth_vel_clipped_fortran_(my_real* vec_move_clippedx, my_real* vec_move_clippedy, my_real* min_pos_Se, my_real *dt);
void get_clipped_ith_vertex_fortran_(long long int *k, Point2D *pt);
void get_clipped_edges_ith_vertex_fortran_(long long int *k, long long int *signed_eR, long long int *signed_eL);
void update_clipped_fortran_(const my_real* vec_move_clippedy, const my_real* vec_move_clippedz, const my_real* dt, \
                            my_real *minimal_length, my_real *maximal_length, my_real *minimal_angle);

#endif