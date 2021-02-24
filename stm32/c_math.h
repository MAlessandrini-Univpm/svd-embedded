#ifndef C_MATH_H
#define C_MATH_H

// Common matrix operations for C code

#include <stdbool.h>


typedef struct {
	float *p;
	int rows, columns, step;
} Mat;


typedef struct {
	int i, j;
	float c, s;
} Rot_info;


void sub_matrix (float *d, const float *s, int rows, int cols, int src_step, bool transpose);
void eye(float*, int n);
void householder(float*, int rows, int columns, float *U, float *V);  // uses AUX1, P, Q
float vector_norm(const float*, int n);
float vector_norm_square(const float*, int n);
void matrix_mul(const Mat *m1, const Mat *m2, Mat *m3);
void matrix_mul_symm(const Mat *m1, Mat *m3);
void matrix_mul_symm_bidiag(const Mat *m1, Mat *m3);
float sign(float n);
float signNo0(float n);
void givens_rotation_left(float c, float s, int i, int k, bool tran, Mat *mat);
void givens_rotation_right(float c, float s, int i, int k, bool tran, Mat *mat);
void golub_kahan_step(float *B, const int n, const int p, const int q, float *U, float *V, const int m);
void eig2(const Mat *m, float *l1, float *l2);
void rotate(float a, float b, float *c, float *s);


#define ONE_SIDED_JACOBI_ROW_MAJOR 0


#endif
