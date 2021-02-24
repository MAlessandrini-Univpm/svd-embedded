#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include "c_math.h"


extern float B[], AUX1[];

static const float eps = 1e-4;
static const float thresh = 1e-6;


int svd_one_sided_jacobi_C(int M, int N, float *U, float *V) {
	// input: B (will be changed), column-major order
	// uses AUX1
	if (V) eye(V, N);  // row-major
	bool exit_flag = false;
	int iterations = 0;
	while (!exit_flag) {
		++iterations;
		exit_flag = true;
		for (int j = N - 1; j >= 1; --j)
			for (int i = j - 1; i >= 0; --i) {
				float alpha = 0, beta = 0, gamm = 0;
#if !ONE_SIDED_JACOBI_ROW_MAJOR
				float *pi = B + M * i, *pj = B + M * j;
				for (int k = 0; k < M; ++k) {
					alpha += *pi * *pi;
					beta += *pj * *pj;
					gamm += *pi++ * *pj++;
				}
#else
				float *pi = B + i, *pj = B + j;
				for (int k = 0; k < M; ++k) {
					alpha += *pi * *pi;
					beta += *pj * *pj;
					gamm += *pi * *pj;
					pi += N;
					pj += N;
				}
#endif
				if (iterations < 50) {
					const float limit = fabsf(gamm) / sqrtf(alpha * beta);
					if (limit > eps) exit_flag = false;
				}
				float c, s;
				if (fabsf(gamm) < thresh) {
					c = 1;
					s = 0;
				}
				else {
					// some computations (square + square root) need to be done in double precision (64 bits)
					// or accuracy does not reach values comparable to other algorithms
					const float tao = (beta - alpha) / (2 * gamm);
					// t can be computed at 32-bit precision, tests show little loss of accuracy
					//  but good speed improvement
					const float t = sign(tao) / (fabsf(tao) + sqrtf(1 + tao * tao));  // t computed at 32-bit precision
					//const double tao64 = tao;
					//const float t = sign(tao) / (fabsf(tao) + (float)sqrt(1 + tao64 * tao64));  // t computed at 64-bit precision
					// tests show that c must instead be computed at 64-bit precision
					//const float c = 1 / sqrtf(1 + t * t);  // c computed at 32-bit precision
					//c = 1 / (float)sqrt(1 + (double)t * (double)t);  // c computed at 64-bit precision
					c = expf(-0.5f * log1pf(t * t));  // new trick by Giorgio! Better than passing to 64 bits.
					s = c * t;
				}
#if !ONE_SIDED_JACOBI_ROW_MAJOR
				// manual Givens rotation of B because it's column-major
				pi = B + M * i; pj = B + M * j;
				for (int k = 0; k < M; ++k) {
					const float t = *pi;
					*pi++ = c * t - s * *pj;
					*pj = s * t + c * *pj;
					++pj;
				}
#else
				Mat mB = {B, M, N, N};
				givens_rotation_right(c, s, i, j, false, &mB);
#endif
				if (V) {
					Mat mV = {V, N, N, N};
					givens_rotation_right(c, s, i, j, false, &mV);
				}
			}
	}
	for (int j = 0; j < N; ++j) {
#if !ONE_SIDED_JACOBI_ROW_MAJOR
		float t = 0, *pj = B + M * j;
		for (int k = 0; k < M; ++k, ++pj) t += *pj * *pj;
#else
		float t = 0, *pj = B + j;
		for (int k = 0; k < M; ++k, pj += N) t += *pj * *pj;
#endif
		AUX1[j] = sqrtf(t);
	}
	if (U) {
		// copy B to U row-major, dividing columns by their singular value
		for (int j = 0; j < N; ++j) {
#if !ONE_SIDED_JACOBI_ROW_MAJOR
			float *pj = B + M * j, *pu = U + j, val = AUX1[j];
			for (int i = 0; i < M; ++i, pu += M) *pu = *pj++ / val;
#else
			float *pj = B + j, *pu = U + j, val = AUX1[j];
			for (int i = 0; i < M; ++i, pu += M, pj += N) *pu = *pj / val;
#endif
		}
	}
	return iterations;
}
