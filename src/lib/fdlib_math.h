#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

__host__ __device__
void
fdlib_math_invert2x2(float matrix[2][2]);

__host__ __device__
void
fdlib_math_invert3x3(float m[][3]);

__host__ __device__
void
fdlib_math_matmul2x2(float A[][2], float B[][2], float C[][2]);

__host__ __device__
void
fdlib_math_matmul3x3(float A[][3], float B[][3], float C[][3]);

__host__ __device__
void
fdlib_math_cross_product(float *A, float *B, float *C);

__host__ __device__
float
fdlib_math_dot_product(float *A, float *B);

__host__ __device__
float
fdlib_math_dist_point2plane(float x0[3], float x1[3], float x2[3], float x3[3]);

__host__ __device__
float
fdlib_math_dist_point2line(float x0, float z0, float p1[2], float p2[2]);

__host__ __device__
void
fdlib_math_bubble_sort(float a[], int index[], int n);

__host__ __device__
void
fdlib_math_bubble_sort_int(int a[], int index[], int n);

__host__ __device__
int
fdlib_math_isPoint2InQuad(float px, float py, const float *vertx, const float *verty);

#endif
