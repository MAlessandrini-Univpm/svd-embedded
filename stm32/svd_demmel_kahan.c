#include <math.h>
#include "c_math.h"
#include <float.h>


extern float B[], AUX1[], AUX2[];  // Z2 cannot be used
static float *lambdas, *mu, lambdaMin, muMin, sigmaLower, sigmaUpper, s1, e1, thresh, maxit;
static const float eps = 1e-6, tol = 1e-3;

static void compute_lambdas_mu(const int n);
static void rot(float a, float b, float *cs, float *sn, float *r);


int svd_demmel_kahan_C(int rows, int columns, float *U, float *V) {
	// input: B (will be changed)
	// uses AUX1, AUX2, P, Q
	householder(B, rows, columns, U, V);  // uses AUX1, P, Q
	const int n = columns;  // n is a size
	maxit = 3 * n * n;
	int iterations = 0;
	while (1) {  // #28
		++iterations;
		bool flag_gk = false;
		// from here: uses AUX1 (lambdas), AUX2 (mu)
		compute_lambdas_mu(n);  // mu size is (n-1)
		for (int j = 0; j < n - 1; ++j)  // #48
			if (fabsf(B[j * n + j + 1] / mu[j]) <= tol) B[j * n + j + 1] = 0;
		for (int j = n - 2; j >= 0; --j)  // #54
			if (fabsf(B[j * n + j + 1] / lambdas[j + 1]) <= tol) B[j * n + j + 1] = 0;
		float tmp = mu[0] / sqrtf(n - 1);  // #60
		for (int i = 1; i < n - 1; ++i)
			if ((mu[i] / sqrtf(n - 1)) < tmp) tmp = mu[i] / sqrtf(n - 1);
		if (B[n * n - n - 1] * B[n * n - n - 1] <= 0.5 * tol * (tmp * tmp - fabsf(B[n * n - 1] * B[n * n - 1])))  // #66
			B[n * n - n - 1] = 0;
		tmp = lambdas[1] / sqrtf(n - 1);  // #70
		for (int i = 2; i < n - 1; ++i) if ((lambdas[i] / sqrtf(n - 1)) < tmp) tmp = lambdas[i] / sqrtf(n - 1);
		if ((B[1] * B[n]) <= 0.5 * tol * (tmp * tmp - fabsf(B[0] * B[0]))) B[1] = 0;  // #76

		int p = -1;
		for (int i = 0; i <= n - 2; ++i) {
			if (fabsf(B[i * n + i + 1]) <= thresh) p = i;
			else break;
			if (i == n - 2) p = n - 1;
		}
		int q = -1;
		for (int i = n - 2; i >= 0; --i) {
			if (fabsf(B[i * n + i + 1]) <= thresh) ++q;
			else break;
			if (i == 0) q = n - 1;
		}

		if (q == n - 1) break;  // end of algorithm  // #108

		for (int i = p + 1; i <= n - q - 3; ++i) {  // #115
			if (B[i * n + i] == 0) {
				for (int k = i + 1; k < n; ++k) {  // #118
					float c, s;
					rotate(B[k * n + k], B[i * n + k], &c, &s);
					Mat m = { B + (p + 1) * n + p + 1, n - p - 1, n - p - 1, n };
					givens_rotation_left(c, s, i, k, true, &m);
				}
			}
			else {
				compute_lambdas_mu(n);  // #126
				if ((n * tol * sigmaLower) <= (eps * sigmaUpper)) {  // #142
					float cs = 1, oldcs = 1, oldsn = 1, sn , r;
					for (int k = p + 1; k <= n - q - 3; ++k) {  // #148
						rot(B[k * n + k] * cs, B[k * n + k + 1], &cs, &sn, &r);
						if (k != p) B[(k - 1) * n + k] = oldsn * r;
						Mat mV = { V, n, n, n };
						if (V) givens_rotation_right(cs, sn, k, k + 1, true, &mV);
						rot(oldcs * r, B[(k + 1) * n + k + 1] * sn, &oldcs, &oldsn, B + k * n + k);
						Mat mU = { U, rows, n, rows };
						if (U) givens_rotation_right(oldcs, oldsn, k, k + 1, true, &mU);
					}
					float h = B[(n - q - 2) * n + n - q - 2] * cs;
					B[(n - q - 3) * n + n - q - 2] = h * oldsn;
					B[(n - q - 2) * n + n - q - 2] = h * oldcs;
				}
				else if (!flag_gk) {
					flag_gk = true;
					golub_kahan_step(B, n, p, q, U, V, rows);  // uses P
				}
			}
		}
	}
	for (int i = 0; i < n; ++i)
		if (B[i * n + i] < 0) {
			B[i * n + i] = -B[i * n + i];
			if (U) for (int j = 0; j < rows; ++j) U[j * rows + i] = -U[j * rows + i];
		}
	return iterations;
}


static void compute_lambdas_mu(const int n) {
	// computes mu, lambda, lambdaMin, muMin, sigmaLower, sigmaUpper, s1, e1, thresh
	// uses AUX1, AUX2
	lambdas = AUX1;
	lambdaMin = lambdas[n - 1] = fabsf(B[n * n - 1]);
	float *p = B + n * (n - 1) - 2;
	for (int j = n - 2; j >= 0; --j, p -= (n + 1)) {
		lambdas[j] = fabsf(p[0]) * lambdas[j + 1] / (lambdas[j + 1] + fabsf(p[1]));
		if (lambdas[j] < lambdaMin) lambdaMin = lambdas[j];
	}
	mu = AUX2;
	mu[0] = fabsf(B[0]);
	p = B + n + 1;
	for (int j = 0; j < n - 1; ++j, p += (n + 1)) {
		mu[j + 1] = fabsf(p[0]) * mu[j] / (mu[j] + fabsf(p[-n]));
		if (j == 0) muMin = mu[j + 1];
		if (mu[j + 1] < muMin) muMin = mu[j + 1];
	}
	++mu;  // from second element to end (n - 1)
	sigmaLower = muMin < lambdaMin ? muMin : lambdaMin;
	p = B + n * (n - 1) - 2;
	s1 = fabsf(B[n * n - 1]);
	e1 = fabsf(p[1]);
	for (int j = 0; j < n - 1; ++j, p -= (n + 1)) {
		if (fabsf(p[0]) > s1) s1 = fabsf(p[0]);
		if (fabsf(p[1]) > e1) e1 = fabsf(p[1]);
	}
	sigmaUpper = s1 > e1 ? s1 : e1;
	thresh = tol * sigmaLower;
	if (maxit * FLT_MIN > thresh) thresh = maxit * FLT_MIN;
}


static void rot(float a, float b, float *cs, float *sn, float *r) {
	if (a == 0) {
		*cs = 0;
		*sn = 1;
		*r = b;
	}
	else if (fabsf(a) > fabsf(b)) {
		const float t = b / a;
		const float tt = expf(0.5f * log1pf(t * t));
		*cs = 1 / tt;
		*sn = t * *cs;
		*r = a * tt;
	}
	else {
		const float t = a / b;
		const float tt = expf(0.5f * log1pf(t * t));
		*sn = 1 / tt;
		*cs = t * *sn;
		*r = b * tt;
	}
}
