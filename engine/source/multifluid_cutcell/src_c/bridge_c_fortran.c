#define _USE_MATH_DEFINES

#include "my_real.h"
#include "bridge_c_fortran.h"
#include "array_double.h"
#include "vector_double.h"
#include "compute_lambdas3D.h"
#include "solver_for_graphblas.h"
#include "update_clipped.h"
# define M_PI           3.14159265358979323846

Polygon2D* grid = NULL;
Polygon2D* clipped = NULL;
Polyhedron3D* clipped3D = NULL;

void launch_grb_(){
    //GrB_init(GrB_NONBLOCKING);
    GrB_init(GrB_BLOCKING);
}

void end_grb_(){
    if(grid) {
        dealloc_Polygon2D(grid); free(grid); grid = NULL;
    }
    if(clipped){
        dealloc_Polygon2D(clipped); free(clipped); clipped = NULL;
    }
    if(clipped3D){
        dealloc_Polyhedron3D(clipped3D); free(clipped3D); clipped3D = NULL;
    }

    GrB_finalize();
}

///  @brief Builds grid given consecutive points with coordinates (x_v, y_v).
///  @details If grid is already built, it only changes the coordinates of the points (since this is always the same structure of edges and faces).
///  @warning It supposes that all grid cells have the same number of points (usually 3 or 4).
void build_grid_from_points_fortran_(const my_real* x_v, const my_real* y_v, const long long int *signed_nb_pts){
    const unsigned long nb_pts = (unsigned long) *signed_nb_pts;
    Point2D p;
    unsigned long i;
    
    if (grid){
        for(i=0; i<nb_pts; i++){
            p = (Point2D){x_v[i], y_v[i]};
            set_ith_elem_vec_pts2D(grid->vertices, i, &p);
        }
    } else {
        grid = polygon_from_consecutive_points(x_v, y_v, nb_pts);
    }
}

/// @brief Take consecutive points to build a polygon.
/// @param x_v array of x-coordinates of size signed_nb_pts
/// @param y_v array of y-coordinates of size signed_nb_pts 
/// @param signed_nb_pts size of the arrays (and number of points).
void build_clipped_from_pts_fortran_(const my_real* x_v, const my_real* y_v, const long long* limits, const long long *nb_polygons){
    unsigned long nb_pts;
    int64_t i;
    long long limit1, limit2;
    const my_real *begin_x, *begin_y;
    Polygon2D* copy_clipped = new_Polygon2D();
    
    if ((*nb_polygons)>0){
        nb_pts = (unsigned long) limits[1];
        clipped = polygon_from_consecutive_points(x_v, y_v, nb_pts);
        for(i=1; i<(*nb_polygons); i++){
            limit1 = limits[i];
            limit2 = limits[i+1]-1;
            nb_pts = (unsigned long)(limit2-limit1+1);
            begin_x = x_v + limit1;
            begin_y = y_v + limit1;
            copy_Polygon2D(clipped, copy_clipped);
            dealloc_Polygon2D(clipped); free(clipped);
            clipped = fuse_polygons(copy_clipped, polygon_from_consecutive_points(begin_x, begin_y, nb_pts));
        }
    } else {
        clipped = new_Polygon2D();
    }
    
    dealloc_Polygon2D(copy_clipped); free(copy_clipped);
}

void compute_lambdas2d_fortran_(const my_real *dt, \
                        my_real *ptr_lambdas_arr,   \
                        my_real *ptr_big_lambda_n,  \
                        my_real *ptr_big_lambda_np1,\
                        Point3D *mean_normal, bool *is_narrowband)
{
    Array_double *lambdas;
    Vector_double *Lambda_n;
    Vector_double *Lambda_np1;

    compute_lambdas2D(grid, clipped3D, *dt, &lambdas, &Lambda_n, &Lambda_np1, mean_normal, is_narrowband);

    memcpy(ptr_lambdas_arr, lambdas->data, lambdas->ncols*lambdas->nrows*sizeof(my_real));

    memcpy(ptr_big_lambda_n, Lambda_n->data, Lambda_n->size*sizeof(my_real));
    memcpy(ptr_big_lambda_np1, Lambda_np1->data, Lambda_np1->size*sizeof(my_real));

    dealloc_arr_double(lambdas); free(lambdas);
    dealloc_vec_double(Lambda_n); free(Lambda_n);
    dealloc_vec_double(Lambda_np1); free(Lambda_np1);
}



void nb_pts_clipped_fortran_(long long int* signed_nb_pts_solid){
    *signed_nb_pts_solid = ((long long) clipped->vertices->size);
}

void nb_edge_clipped_fortran_(long long int* signed_nb_edges_solid){
    unsigned long nb;
    GrB_Matrix_ncols(&nb, *(clipped->edges));
    *signed_nb_edges_solid = ((long long) nb);
}

void compute_normals_clipped_fortran_(my_real* normalVecx, my_real* normalVecy, \
                                        my_real* normalVecEdgex, my_real* normalVecEdgey, \
                                        my_real* min_pos_Se){

    Vector_points2D* normals_pts;
    Vector_points2D* normals_edges;
    GrB_Index nb_pts, nb_edges;
    uint64_t i;

    GrB_Matrix_nrows(&nb_pts, *(clipped->edges));
    GrB_Matrix_ncols(&nb_edges, *(clipped->edges));
    normals_pts = alloc_with_capacity_vec_pts2D(nb_pts);
    normals_pts->size = nb_pts;
    normals_edges = alloc_with_capacity_vec_pts2D(nb_edges);
    normals_edges->size = nb_edges;
    compute_all_normals2D(clipped, normals_pts, normals_edges, min_pos_Se);


    for(i=0; i<normals_pts->size; i++){
        normalVecx[i] = normals_pts->points[i].x;
        normalVecy[i] = normals_pts->points[i].y;
    }
    for(i=0; i<normals_edges->size; i++){
        normalVecEdgex[i] = normals_edges->points[i].x;
        normalVecEdgey[i] = normals_edges->points[i].y;
    }

    dealloc_vec_pts2D(normals_pts); free(normals_pts);
    dealloc_vec_pts2D(normals_edges); free(normals_edges);
}

void smooth_vel_clipped_fortran_(my_real* vec_move_clippedx, my_real* vec_move_clippedy, my_real* min_pos_Se, my_real *dt){
    my_real eps;
    GrB_Matrix smoothing_op, id;
    unsigned long nb_pts = clipped->vertices->size;
    GrB_Index ind_vector[2] = {1, nb_pts};
    
    eps = 0.001 / ((*min_pos_Se)*(*min_pos_Se));
    GrB_Matrix_new(&smoothing_op, GrB_FP64, nb_pts, nb_pts);
    GrB_Matrix_new(&id, GrB_FP64, nb_pts, nb_pts);

    //smoothing_op = Id + eps*dt*(edges*edges')
    GrB_mxm(smoothing_op, GrB_NULL, GrB_NULL, GrB_NULL, *(clipped->edges), *(clipped->edges), GrB_DESC_T1);
    GrB_assign(id, GrB_NULL, GrB_NULL, 1.0, ind_vector, GxB_RANGE, ind_vector, GxB_RANGE, GrB_NULL);
    GrB_apply(smoothing_op, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, eps*(*dt), smoothing_op, GrB_NULL);
    GrB_eWiseAdd(smoothing_op, GrB_NULL, GrB_NULL, GrB_PLUS_FP64, id, smoothing_op, GrB_NULL);

    solver_for_graphblas_FP_sym(smoothing_op, vec_move_clippedx);
    solver_for_graphblas_FP_sym(smoothing_op, vec_move_clippedy);

    GrB_free(&smoothing_op);
    GrB_free(&id);
}

void get_clipped_ith_vertex_fortran_(long long *signed_k, Point2D *pt){
    uint64_t k = (uint64_t)((*signed_k) - 1);
    *pt = *get_ith_elem_vec_pts2D(clipped->vertices, k); 
}

void get_clipped_edges_ith_vertex_fortran_(long long *k_signed, long long *signed_eR, long long *signed_eL){
    uint64_t eR, eL;
    GrB_Matrix e_k;
    GrB_Vector nz_e_k, extr_vals_e_k, I_vec_e_k;
    GrB_Index size_nz_e_k;
    GrB_Index nb_pts;
    GrB_Info infogrb;
    GrB_Index k = (GrB_Index)(*k_signed - 1);

    nb_pts = clipped->vertices->size;

    infogrb = GrB_Matrix_new(&e_k, GrB_INT8, 1, nb_pts);
    infogrb = GrB_Vector_new(&nz_e_k, GrB_UINT64, 2);
    infogrb = GrB_Vector_new(&extr_vals_e_k, GrB_INT8, 2);
    infogrb = GrB_Vector_new(&I_vec_e_k, GrB_UINT64, 2);
    GrB_extract(e_k, GrB_NULL, GrB_NULL, *(clipped->edges), &k, 1, GrB_ALL, 1, GrB_NULL); //Get indices of edges connected with point k
    GxB_Matrix_extractTuples_Vector(I_vec_e_k, nz_e_k, extr_vals_e_k, e_k, GrB_NULL);
    GrB_Vector_size(&size_nz_e_k, nz_e_k);
    if (size_nz_e_k == 2){ //Exactly two edges connected to the point k
        GrB_Vector_extractElement(&eR, nz_e_k, 0);
        GrB_Vector_extractElement(&eL, nz_e_k, 1);
        *signed_eR = (long long)eR;
        *signed_eL = (long long)eL;
    } else {                   
        *signed_eR = -1;       
        *signed_eL = -1;       
    }                          
    GrB_free(&e_k);            
    GrB_free(&nz_e_k);         
    GrB_free(&extr_vals_e_k);  
    GrB_free(&I_vec_e_k);      
}                              
                               
void update_clipped_fortran_(const my_real* vec_move_clippedy, const my_real* vec_move_clippedz, const my_real* dt, \
                            my_real *minimal_length, my_real *maximal_length, my_real *minimal_angle){

    update_solid(&clipped, &clipped3D, vec_move_clippedy, vec_move_clippedz, *dt, \
                    *minimal_length, *maximal_length, *minimal_angle);
}                              
                               
                               
                               
                               
                               