#include <math.h>
#include "c_math.h"


extern float B[];  // Z2 cannot be used
extern int I1[];
static const float eps = 1e-4;
//extern uint32_t precisionTimer;


int svd_golub_reinsch_C(int rows, int columns, float *U, float *V) {
	// input: B (will be changed)
	// uses AUX1, P, Q, I1
	//precisionTimer = 0;
	householder(B, rows, columns, U, V);  // uses AUX1, P, Q
	const int n = columns;  // n is a size
	int *nulleigindx = I1;
	int nulleigindx_sz = 0;
	int iterations = 0;
	while(1) {
		++iterations;
		// Set superdiagonal elements to zero if less than precision........
		for (int i = 0; i <= n - 2; ++i)
			if (fabsf(B[i * n + i + 1]) < eps * (fabsf(B[i * n + i]) + fabsf(B[(i + 1) * n + i + 1])))
				B[i * n + i + 1] = 0;
		// p & q evaluation............................................
		int p = -1;
		for (int i = 0; i <= n - 2; ++i) {
			if (fabsf(B[i * n + i + 1]) < eps) p = i;
			else break;
			if (i == n - 2) p = n - 1;
		}
		int q = -1;
		for (int i = n - 2; i >= 0; --i) {
			if (fabsf(B[i * n + i + 1]) < eps) ++q;
			else break;
			if (i == 0) q = n - 1;
		}
		// Stop if null superdiagonal is found..................
		if (q == n - 1)	break;
		// Givens rotations in the case some eigenvalues are 0............
		bool GK_step = true;
		for (int i = p + 1; i <= n - q - 3; ++i) {
			// find if i not in nulleigindx
			bool newzeroindx = true;
			for (int j = 0; j < nulleigindx_sz; ++j) if (nulleigindx[j] == i) newzeroindx = false;
			if (fabsf(B[i * n + i]) < eps && newzeroindx) {
				nulleigindx[nulleigindx_sz++] = i;
				GK_step = false;
				for (int k = i + 1; k <= n - q - 2; ++k) {
					float c, s;
					rotate(B[k * n + k], B[i * n + k], &c, &s);
					Mat m = { B, n, n, n };
					givens_rotation_left(c, s, i, k, true, &m);
				}
			}
		}
		if (GK_step) {
			//systick_timer_start();
			golub_kahan_step(B, n, p, q, U, V, rows);  // uses P
			//precisionTimer += systick_timer_stop();
		}
	}
	for (int i = 0; i < n; ++i)
		if (B[i * n + i] < 0) {
			B[i * n + i] = -B[i * n + i];
			if (U) for (int j = 0; j < rows; ++j) U[j * rows + i] = -U[j * rows + i];
		}
	return iterations;
}
