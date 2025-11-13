#include "compute_lambdas3D.h"
#include <math.h>
#include "update_clipped.h"

//Clips `p` using a plane defined by `n_face` and `pt_face`.
//`mark_edge` is used as a tag for retrieving the edge or face defining the cutting plane.
//`sign_taken` changes the orientation of `n_face`.
static void build_clipped_in_partial(Polyhedron3D* p, Point3D *n_face, Point3D *pt_face, long int mark_edge, int8_t sign_taken){
    Vector_points3D *pts_copy = alloc_empty_vec_pts3D();
    Vector_int *status_face = alloc_empty_vec_int();
    GrB_Matrix *edges_in, *faces_in, *volumes_in;
    
    edges_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    faces_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));
    volumes_in = (GrB_Matrix*)malloc(sizeof(GrB_Matrix));

    cut_edges3D(p, n_face, pt_face, sign_taken, pts_copy, edges_in);

    close_cells(edges_in, p->faces, NULL, -1, faces_in);

    copy_vec_int(p->status_face, status_face);
    close_cells(faces_in, p->volumes, status_face, mark_edge, volumes_in);

    //Copy new polyhedron in p.
    copy_vec_pts3D(pts_copy, p->vertices);
    dealloc_vec_pts3D(pts_copy); free(pts_copy);
    copy_vec_int(status_face, p->status_face);
    dealloc_vec_int(status_face); free(status_face);
    GrB_Matrix_dup(p->edges, *edges_in);
    GrB_Matrix_dup(p->faces, *faces_in);
    GrB_Matrix_dup(p->volumes, *volumes_in);
    GrB_free(edges_in); free(edges_in);
    GrB_free(faces_in);free(faces_in);
    GrB_free(volumes_in);free(volumes_in);
}

//Returns the first point of `i`th face in p.
static Point3D* first_point_of_ith_face(const Polyhedron3D *p, GrB_Index i){
    Point3D* pt;
    GrB_Index nb_faces, e0, pt_index;
    GrB_Info infogrb;
    GrB_Vector whole_vec, nz, extr_vals;
    GrB_Index nrows;

    GrB_Matrix_ncols(&nb_faces, *(p->faces));
    if (i>=nb_faces){
        pt = (Point3D*)malloc(sizeof(Point3D));
        pt->x = 0./0.;
        pt->y = 0./0.;
        pt->t = 0./0.;
    } else {
        infogrb = GrB_Matrix_nrows(&nrows, *(p->faces));
        infogrb = GrB_Vector_new(&whole_vec, GrB_INT8, nrows);
        infogrb = GrB_Vector_new(&nz, GrB_UINT64, nrows);
        infogrb = GrB_Vector_new(&extr_vals, GrB_INT8, nrows);
        infogrb = GrB_extract(whole_vec, GrB_NULL, GrB_NULL, *(p->faces), GrB_ALL, 1, i, GrB_NULL); //Get indices of edges composing face i
        infogrb = GxB_Vector_extractTuples_Vector(nz, extr_vals, whole_vec, GrB_NULL);
        infogrb = GrB_Vector_extractElement(&e0, nz, 0);//Get index of first edge composing face i

        infogrb = GrB_Matrix_nrows(&nrows, *(p->edges));
        infogrb = GrB_Vector_resize(whole_vec, nrows);
        infogrb = GrB_Vector_resize(nz, nrows);
        infogrb = GrB_Vector_resize(extr_vals, nrows);
        infogrb = GrB_extract(whole_vec, GrB_NULL, GrB_NULL, *(p->edges), GrB_ALL, 1, e0, GrB_NULL); //Get indices of points of first edge composing face i
        infogrb = GxB_Vector_extractTuples_Vector(nz, extr_vals, whole_vec, GrB_NULL);
        infogrb = GrB_Vector_extractElement(&pt_index, nz, 0);//Get index of first point of first edge composing face i

        pt = get_ith_elem_vec_pts3D(p->vertices, pt_index);
    }

    GrB_free(&whole_vec);
    GrB_free(&nz);
    GrB_free(&extr_vals);
    return pt;
}

/// @brief Clips `clipped` using `clipper`.
/// @details We suppose that `clipper` is convex and that `clipper.status_face < 0` when the corresponding face is at t==t^n or t==t^{n+1}.
/// @param clipper [IN] a convex polygon
/// @param clipped [IN] the clipped polygon
/// @return The result of the clipping (a polyhedron)
Polyhedron3D* clip3D(const Polyhedron3D *clipper, const Polyhedron3D *clipped){
    Polyhedron3D* clipped_in = new_Polyhedron3D();
    Vector_points3D* normal_vectors = points3D_from_matrix(surfaces_poly3D(clipper));
    GrB_Index i;//, size_nv;
    Point3D *n_face, *pt_face;
    long int *stf_i;
    int8_t vol_i;

    copy_Polyhedron3D(clipped, clipped_in);
    //for i in eachindex(normal_vectors)
    for (i = 0; i<normal_vectors->size; i++){
        if(!clipped_in)
            break;
        stf_i = get_ith_elem_vec_int(clipper->status_face, i);
        if (*stf_i > 2){
            n_face = get_ith_elem_vec_pts3D(normal_vectors, i);
            pt_face = first_point_of_ith_face(clipper, i);
            GrB_Matrix_extractElement(&vol_i, *(clipper->volumes), i, 0);
            build_clipped_in_partial(clipped_in, n_face, pt_face, *stf_i, -vol_i);
            clean_Polyhedron3D(clipped_in, &clipped_in);
        }
    }

    dealloc_vec_pts3D(normal_vectors); free(normal_vectors);

    return clipped_in;
}

/// @brief Compute the effective area on each face of grid.
/// @details The effective area is the area of each face of `grid` minus the area intersected with `initial_p`.
///          The result consists in an ordered list of the faces of `grid` and 
///          `lambdas`, the effective area in a list ordered in the same order.
///          In `lambdas`, we store the effective area and the "complementarity", which is the area of the face minus the effective area.
/// @param grid [IN] Clipper cell
/// @param initial_p [IN] Clipped polygon
/// @param lambdas [OUT] effective areas (allocated inside the function)
/// @param mean_normal [OUT] normal vector of the surface of `initial_p` clipped inside `grid`.
/// @param is_narrowband [OUT] true if the intersection of `grid` and `initial_p` is not empty, false otherwise.
void compute_lambdas2D_time(const Polyhedron3D* grid, const Polyhedron3D *initial_p, \
                            Vector_points3D **lambdas, Point3D *mean_normal, bool *is_narrowband){
    Polyhedron3D *p = clip3D(grid, initial_p);
    GrB_Index nb_edge, nb_cols_vol, nb_cols_fac, i, j;//, e;
    Point3D *nvpi, *lam;
    int8_t pvij_int;
    long int *psfi;
    my_real pvij;
    Vector_points3D *norm_vec_poly;

    GrB_Matrix_ncols(&nb_edge, *(grid->faces)); //actually number of edges + 2 faces at times tn and tn+dt

    *lambdas = alloc_with_capacity_vec_pts3D(nb_edge);
    *mean_normal = (Point3D){0.0, 0.0, 0.0};
    for (j=0; j<nb_edge; j++){
        set_ith_elem_vec_pts3D(*lambdas, j, mean_normal);
    }

    *is_narrowband = false;
    
    if(p){
        norm_vec_poly = points3D_from_matrix(surfaces_poly3D(p));
    
        GrB_Matrix_ncols(&nb_cols_vol, *(p->volumes)); 
        GrB_Matrix_ncols(&nb_cols_fac, *(p->faces)); 
        for(i=0; i<nb_cols_fac; i++){
            psfi = get_ith_elem_vec_int(p->status_face, i);
            nvpi = get_ith_elem_vec_pts3D(norm_vec_poly, i);
            for (j=0; j<nb_cols_vol; j++){
                pvij_int = 0; //if volumes[i,j] does not exist, it won't change the value of pvij_int.
                GrB_Matrix_extractElement(&pvij_int, *(p->volumes), i, j);
                if (pvij_int != 0){
                    if (*psfi<1){
                        pvij = (my_real) pvij_int;
                        inplace_xpay_points3D(mean_normal, pvij, nvpi);
                        *is_narrowband = true;
                    }
                    else {
                        lam = get_ith_elem_vec_pts3D(*lambdas, *psfi-1);
                        inplace_axpy_points3D(lam, 1.0, nvpi);
                    }
                    //for (e = 1; e<=nb_edge; e++){
                    //    if (*psfi == ((long int) e)){
                    //        lam = get_ith_elem_vec_pts3D(lambdas, e-1);
                    //        inplace_axpy_points3D(lam, 1.0, nvpi);
                    //    }
                    //}
                }
            }
        }
    }

    if(p){
        dealloc_vec_pts3D(norm_vec_poly); free(norm_vec_poly);
        dealloc_Polyhedron3D(p); free(p);
    }
}

/// @brief Compute the effective areas in `grid` when occupied by `clipped` moving with vectors `vec_move_solid` for time `dt`. 
/// @details vec_move_solid should be an array of array of size #columns(clipped.faces), 
///           and each member i of this vector should be of size mini_clipped->vertices->size, 
///           where mini_clipped = extract_ith_face(clipped, i)
/// @param grid [IN] non-moving polygon, clipper.
/// @param clipped3D [IN] moving polygon.
/// @param dt [IN] time-step
/// @param nb_pts [IN] Number of points in 2D polygon used to build clipped3D
/// @param id_cell [IN] Cell identificator of grid
/// @param lambdas_arr [OUT] Array of effective area for each edge of `grid`. Allocated inside the function.
/// @param big_lambda_n [OUT] First cell: effective area at time t^n. Second cell: area of grid - effective area. Allocated inside the function.
/// @param big_lambda_np1 [OUT] First cell: effective area at time t^n+`dt`. Second cell: area of grid - effective area. Allocated inside the function.
/// @param is_narrowband [OUT] True if the intersection of `grid` and `clipped` is not empty at time t^n or t^n+`dt`, false otherwise.
void compute_lambdas2D(const Polygon2D* grid, const Polyhedron3D *clipped3D, const my_real dt, \
                        Array_double **lambdas_arr, Vector_double** big_lambda_n, Vector_double** big_lambda_np1, \
                        Point3D *mean_normal, bool *is_narrowband){
    const unsigned int nb_regions = 2;
    my_real *val = (my_real*)malloc(sizeof(my_real));
    my_real nm;
    Point3D *area, *occupied;
    Vector_points3D *occupied_area;
    //Polygon2D *mini_clipped;
    Polyhedron3D *cell3D = NULL;
    //long *sfj;
    Point3D *pt3D, local_mean_normal, *local_l;
    bool local_narrowband;
    //Vector_points2D *vec_move_grid;
    my_real *vec_move_gridx = NULL, *vec_move_gridy = NULL;
    Vector_points3D *surfaces = NULL;
    GrB_Index nb_edge, i, j, k;
    //GrB_Index nb_clipped_faces;
    Vector_points3D *lambdas3D; //TODO : Change this when nb_regions>2
    Array_points3D *local_lambdas;
    
    GrB_Matrix_ncols(&nb_edge, *(grid->edges));
    *lambdas_arr = alloc_with_capacity_arr_double(nb_edge, nb_regions); //All set to 0
    local_lambdas = alloc_with_capacity_arr_pts3D(nb_edge + 2, nb_regions-1); //All set to 0
    occupied_area = alloc_with_capacity_vec_pts3D(nb_edge + 2);

    *big_lambda_n = alloc_with_capacity_vec_double(nb_regions);
    *big_lambda_np1 = alloc_with_capacity_vec_double(nb_regions);
    *val = 0.;
    for (i=0; i<nb_regions; i++){
        set_ith_elem_vec_double(*big_lambda_n, i, val);
        set_ith_elem_vec_double(*big_lambda_np1, i, val);
    }
    
    *mean_normal = (Point3D){0.,0.,0.};
    for (i=0; i<nb_edge + 2; i++){
        set_ith_elem_vec_pts3D(occupied_area, i, mean_normal);
    }

    *is_narrowband = false;
    vec_move_gridx = calloc(grid->vertices->size, sizeof(my_real));
    vec_move_gridy = calloc(grid->vertices->size, sizeof(my_real));
    cell3D = build_space2D_time_cell(grid, vec_move_gridx, vec_move_gridy, grid->vertices->size, dt, false, NULL);
    surfaces = points3D_from_matrix(surfaces_poly3D(cell3D));

    if ((clipped3D) && (clipped3D->vertices->size>2)){
        k = 0; //should be k = clipped->status_edge[i], or another variable to indicate what region covers face nb i.
        
        compute_lambdas2D_time(cell3D, clipped3D, &lambdas3D, &local_mean_normal, &local_narrowband);

        mean_normal->x += local_mean_normal.x;
        mean_normal->y += local_mean_normal.y;
        mean_normal->t += local_mean_normal.t;
        *is_narrowband |= local_narrowband;
        
        //λe += λ[2:end]
        for (j=0; j<lambdas3D->size; j++){
            pt3D = get_ith_elem_vec_pts3D(lambdas3D, j); 
            local_l = get_ijth_elem_arr_pts3D(local_lambdas, j, k);
            local_l->x += pt3D->x;
            local_l->y += pt3D->y;
            local_l->t += pt3D->t;

            local_l = get_ith_elem_vec_pts3D(occupied_area, j);
            local_l->x += pt3D->x;
            local_l->y += pt3D->y;
            local_l->t += pt3D->t;
        }
        

        i = 0;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied));

        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }
        

        i = 1;
        area = get_ith_elem_vec_pts3D(surfaces, i);
        occupied = get_ith_elem_vec_pts3D(occupied_area, i);
        *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            area = get_ith_elem_vec_pts3D(surfaces, i);
            occupied = get_ith_elem_vec_pts3D(occupied_area, i);
            *val = fmax(0., norm_pt3D(*area) - norm_pt3D(*occupied));
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);
            for(k = 1; k<nb_regions; k++){
                nm = norm_pt3D(*get_ijth_elem_arr_pts3D(local_lambdas, i, k-1));
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
    } else {
        i = 0;
        *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_n, 0, val);
        nm = 0.;
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_n, k, &nm);
        }

        i = 1;
        *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i));
        set_ith_elem_vec_double(*big_lambda_np1, 0, val);
        for(k = 1; k<nb_regions; k++){
            set_ith_elem_vec_double(*big_lambda_np1, k, &nm);
        }

        for (i=2; i<nb_edge + 2; i++){
            *val = norm_pt3D(*get_ith_elem_vec_pts3D(surfaces, i));
            set_ijth_elem_arr_double(*lambdas_arr, i-2, 0, val);

            for(k = 1; k<nb_regions; k++){
                set_ijth_elem_arr_double(*lambdas_arr, i-2, k, &nm); 
            }
        }
    }

    dealloc_Polyhedron3D(cell3D); free(cell3D);
    dealloc_vec_pts3D(surfaces); free(surfaces);
    dealloc_vec_pts3D(occupied_area); free(occupied_area);
    dealloc_vec_pts3D(lambdas3D); free(lambdas3D);
    dealloc_arr_pts3D(local_lambdas); free(local_lambdas);
    if (vec_move_gridx) free(vec_move_gridx);
    if (vec_move_gridy) free(vec_move_gridy);
    if (val) free(val);
}
