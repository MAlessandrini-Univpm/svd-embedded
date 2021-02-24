#include <math.h>
#include "c_math.h"
#include <stddef.h>
#include <assert.h>


extern float B[], Q[], AUX1[];

static const float eps = 1e-6;


int svd_jacobi_rotation_C(int rows, int columns, float *U, float *boolRight) {
	// input: B (will be changed)
	// uses AUX1, P, Q
	const int n = columns;  // n is a size
	Mat m1 = { B, n, n, n }, mQ = { Q, n, n, n }, mU = { U, rows, columns, rows };
	for (int i = 0; i < n * n; ++i) assert(!isnan(B[i]) && !isinf(B[i]));
	if (boolRight) {
		householder(B, rows, columns, NULL, U);  // uses AUX1, P, Q
		mU.rows = columns;
		mU.step = columns;
	}
	else {
		householder(B, rows, columns, U, NULL);  // uses AUX1, P, Q
	}
	for (int i = 0; i < n * n; ++i) assert(!isnan(B[i]) && !isinf(B[i]));
	matrix_mul_symm_bidiag(&m1, &mQ);
	// from here we use Q instead of B
	float NN = 0;
	for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) NN += Q[i * n + j] * Q[i * n + j];
	assert(NN > 0);
	int iterations = 0;
	while (1) {
		++iterations;
		float a = 0;
		for (int i = 0; i <= n - 2; ++i)
			for (int j = i + 1; j <= n - 1; ++j) {
				a += Q[i * n + j] * Q[i * n + j] + Q[j * n + i] * Q[j * n + i];
				float c = 1, s = 0;
				if (fabsf(Q[i * n + j]) > eps) {
					float tau = (Q[j * n + j] - Q[i * n + i]) / (2 * Q[i * n + j]);
					float t = sign(tau) / (fabsf(tau) + sqrtf(1 + tau * tau));
					c = expf(-0.5f * log1pf(t * t));
					s = c * t;
				}
				givens_rotation_left(c, s, i, j, true, &mQ);
				givens_rotation_right(c, s, i, j, false, &mQ);
				if (U) givens_rotation_right(c, s, i, j, false, &mU);
			}
		if (a <= eps * eps * NN) break;
	}
	for (int i = 0; i < n; ++i) AUX1[i] = sqrtf(Q[i * n + i]);
	return iterations;
}
