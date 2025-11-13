#ifndef POINT_H
#define POINT_H

#include "my_real.h"

typedef struct 
{
    my_real x, y;
} Point2D;

typedef struct 
{
    my_real x, y, t;
} Point3D;

typedef struct 
{
    my_real x, y, z, t;
} Point4D;

my_real scalProd2D(const Point2D v1, const Point2D v2);
my_real compute_distance2D(const Point2D v, const Point2D n, const Point2D pt);
Point2D find_intersection2D(const Point2D v1, const my_real d1, const Point2D v2, const my_real d2);
void sum_points2D(void* out_void, const void* v1_void, const void* v2_void);
void scale_points2D(void* out_void, const void* s_void, const void* v_void);
void inplace_axpy_points2D(Point2D* v1, const my_real a, const Point2D* v2); //v1 = a*v1+v2;
void inplace_xpay_points2D(Point2D* v1, const my_real a, const Point2D* v2); //v1 = v1+a*v2;
void print_pt2D(const Point2D *p);
my_real norm_pt2D(const Point2D p);

my_real scalProd3D(const Point3D v1, const Point3D v2);
my_real compute_distance3D(const Point3D v, const Point3D n, const Point3D pt);
Point3D find_intersection3D(const Point3D v1, const my_real d1, const Point3D v2, const my_real d2);
void sum_points3D(void* out_void, const void* v1_void, const void* v2_void);
void scale_points3D(void* out_void, const void* s_void, const void* v_void);
void inplace_axpy_points3D(Point3D* v1, const my_real a, const Point3D* v2); //v1 = a*v1+v2;
void inplace_xpay_points3D(Point3D* v1, const my_real a, const Point3D* v2); //v1 = v1+a*v2;
void print_pt3D(const Point3D *p);
my_real norm_pt3D(const Point3D p);

my_real scalProd4D(const Point4D v1, const Point4D v2);
my_real compute_distance4D(const Point4D v, const Point4D n, const Point4D pt);
Point4D find_intersection4D(const Point4D v1, const my_real d1, const Point4D v2, const my_real d2);
void sum_points4D(void* out_void, const void* v1_void, const void* v2_void);
void scale_points4D(void* out_void, const void* s_void, const void* v_void);
void inplace_axpy_points4D(Point4D* v1, const my_real a, const Point4D* v2); //v1 = a*v1+v2;
void inplace_xpay_points4D(Point4D* v1, const my_real a, const Point4D* v2); //v1 = v1+a*v2;
void print_pt4D(const Point4D *p);
my_real norm_pt4D(const Point4D p);

#endif