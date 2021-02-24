#ifndef ALG_H
#define ALG_H

#include <stdint.h>

extern float population[];
float population_mean(unsigned count);
float population_stddev(unsigned count, float max_deviation);

float mean_error(float const *series1a, float const *series1b, unsigned count1,
	float const *series2a, float const *series2b, unsigned count2);  // series2a/b can be NULL

int compare_float_desc(const void *a, const void *b);  // for qsort function

void timer_start(void);
float timer_stop(void);
// precision time measurement using SysTick
void systick_timer_start(void);
uint32_t systick_timer_stop(void);

// SVD functions return the number of iterations performed (where applicable)
// C SVD functions take U, V optionals.
//  U is rows*rows but the result matrix can be used as the rows*cols portion
//  V is columns*columns
int svd_golub_reinsch_C(int, int, float *U, float *V);  // B ← diagonal matrix RxC with unsorted singular values (modifies B, P, Q, AUX1, I1).
int svd_one_sided_jacobi_C(int, int, float *U, float *V);  // AUX1 ← unsorted singular values (modifies B)
int svd_one_sided_jacobi_asm(int, int);  // AUX1 ← unsorted singular values (modifies B)
int svd_jacobi_rotation_C(int, int, float *U, float *boolRight);  // AUX1 ← unsorted singular values (modifies B, P, Q)
int svd_demmel_kahan_C(int, int, float *U, float *V);  // B ← diagonal matrix RxC with unsorted singular values (modifies B, P, Q, AUX1, AUX2)
int svd_divide_conquer_C(int, int, float *U, float *boolRight);  // AUX1 ← unsorted singular values (modifies B, P, Q, U, V, AUX2, Z1, Z3, I1, I2)


#endif
