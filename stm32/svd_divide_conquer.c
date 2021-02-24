#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "c_math.h"


extern float B[], AUX1[], AUX2[], P[], Q[], Z1[], U[], UU[], V[], Z3[];  // Z2 cannot be used
extern int I1[], I2[];
extern Rot_info Rot[];

static const float threshold = 1e-6;

static void dandcOpt(Mat *mT, Mat *mQ, float *mdlam);
static int deflate(const float*, float*, const int, float*, float*, int*, Rot_info*);
static void gu_eisenstat(float*, float*, const int, float*, float*, const int);
static void sort_indexes(int*, const float*, int);
static float zerodandc(/*const*/ float*, const float*, const int, const int, float*);


int svd_divide_conquer_C(int rows, int columns, float *Uout, float *boolRight) {
	// input: B (will be changed)
	// output:
	//  Q: eigenvectors
	//  AUX1: eigenvalues
	// uses:
	//  P: matrix T (tridiagonal symmetric, input to divide and conquer)
	if (boolRight) householder(B, rows, columns, NULL, Uout ?  UU : NULL);
	else householder(B, rows, columns, Uout ? UU : NULL, NULL);
	const int n = columns;
	Mat m1 = { B, n, n, n }, mT = { P, n, n, n };
	matrix_mul_symm_bidiag(&m1, &mT);
	Mat mQ = { Q, n, n, n };
	dandcOpt(&mT, &mQ, AUX1);  // compute eigenvalues/eigenvectors recursively
	for (int i = 0; i < n; ++i) AUX1[i] = sqrtf(AUX1[i]);
	if (Uout) {
		if (boolRight) {
			Mat mVh = { UU, columns, columns, columns };
			Mat mU  = { Uout, columns, columns, columns };
			matrix_mul(&mVh, &mQ, &mU);
		}
		else {
			Mat mUh = { UU, rows, columns, rows };
			Mat mU  = { Uout, rows, columns, rows };
			matrix_mul(&mUh, &mQ, &mU);
		}
	}
	return 0;
}


static void dandcOpt(Mat *mT, Mat *mQ, float *mdlam) {
	// input:
	//  mT (will be changed)
	// output:
	//  mQ, mdlam
	// uses: AUX2, Z1, Z3, I1, Rot, U
	const int n = mT->rows;
	if (n == 1) {
		*mQ->p = 1;
		*mdlam = *mT->p;
		return;
	}
	const int x = n / 2 + n % 2;  // to have the same results as matlab
	Mat mT1 = { mT->p, x, x, mT->step }, mT2 = { mT->p + x * mT->step + x, n - x, n - x, mT->step };
	float rho = *(mT->p + (x - 1) * mT->step + x), mult = 1;
	if (rho > 0) {
		*(mT1.p + (x - 1) * mT1.step + x - 1) -= rho;
		*mT2.p -= rho;
	}
	else {
		*(mT1.p + (x - 1) * mT1.step + x - 1) += rho;
		*mT2.p += rho;
		rho = -rho;
		mult = -1;
	}
	// recursive call
	Mat mQ1 = { mQ->p, x, x, mQ->step }, mQ2 = { mQ->p + x * mQ->step + x, n - x, n - x, mQ->step };
	dandcOpt(&mT1, &mQ1, mdlam);
	dandcOpt(&mT2, &mQ2, mdlam + x);
	// after the recursions have finished, we can reuse mT->p

	float *v = AUX2;
	for (int i = 0; i < x; ++i) v[i] = mult * sqrtf(rho) * *(mQ1.p + (x - 1) * mQ1.step + i);
	for (int i = 0; i < (n - x); ++i) v[x + i] = sqrtf(rho) * *(mQ2.p + i);

	// deflation
	float *v1 = Z1, *d1 = Z3;
	// d1 will contain original d1 + original d2 (elements of d reorganized after deflation)
	// I1 = null_indexes
	// v will be changed, too
	const int sv1 = deflate(mdlam, v, n, v1, d1, I1, Rot);
	const int sv2 = n - sv1;
	//if (sv2) printf("DEFLATION: %d -> %d\n", n, sv1);
	// v is no more used, AUX2 is free
	// mdlam is no more used, will be assigned to solution of gu_eisenstat
	float *Qgu = U;  // n * n (Q1 and Q1_merged in matlab code, solution of gu-eisenstat)
	// let gu_eisenstat fill the top-left portion of Q_gu (undeflated part)
	if (sv1) gu_eisenstat(d1, v1, sv1, Qgu, mdlam, n);
	// d1 (first part of d) and v1 will not be used anymore. d2 will stay intact
	if (sv2 > 0) {
		// augment Qgu and mdlam with the deflated-out parts, obtaining Q1_merged and dlam_merged
		for (int i = 0; i < sv1; ++i) memset(Qgu + i * n + sv1, 0, sv2 * sizeof(float));
		for (int i = sv1; i < n; ++i) memset(Qgu + i * n, 0, n * sizeof(float));
		for (int i = sv1; i < n; ++i) Qgu[i * n + i] = 1;
		for (int i = sv1; i < n; ++i) mdlam[i] = d1[i];
		float *TMP = mT->p;  // reuse mT
		// permute back Qgu (to TMP) reversing deflation
		const float *p1 = Qgu, *p2 = Qgu + sv1 * n;
		const int *nulls = I1;
		for (int i = 0; i < n; ++i) {
			// is i a null index from deflation?
			if (*nulls >= 0 && *nulls == i) {  // yes
				++nulls;
				memcpy(TMP + i * n, p2, n * sizeof(float));
				p2 += n;
			}
			else {
				memcpy(TMP + i * n, p1, n * sizeof(float));
				p1 += n;
			}
		}
		// apply inverse Givens rotation
		for (Rot_info *rot = Rot; rot->i != rot->j; ++rot) {
			const float c = rot->c, s = -rot->s;
			float *r1 = TMP + n * rot->i, *r2 = TMP + n * rot->j;
			for (int k = 0; k < n; ++k, ++r1, ++r2) *r1 = c * *r1 + s * *r2;
			r1 = TMP + n * rot->i; r2 = TMP + n * rot->j;
			for (int k = 0; k < n; ++k, ++r1, ++r2) *r2 = -s * *r1 + c * *r2;
		}
		// final sorting is not necessary
		Qgu = TMP;  // so the pointer is the same even in no deflation happened
	}
	// Now we have to multiply [ mQ1 0 ; 0 mQ2 ] * Qgu and store it in mQ (same memory as mQ1, mQ2).
	// To do the multiplication in-place over mQ we have to save a temporary copy of every (half) row of mQ
	// first half (upper)
	float *r3 = mQ->p;  // also destination
	for (int r = 0; r < x; ++r, r3 += mQ->step) {
		memcpy(AUX2, r3, x * sizeof(float));
		memset(r3, 0, n * sizeof(float));
		float *p3 = r3;
		for (int c = 0; c < n; ++c, ++p3) {
			float *p1 = AUX2, *p2 = Qgu + c;
			for (int k = 0; k < x; ++k, ++p1, p2 += n)
				*p3 += *p1 * *p2;
		}
	}
	// second half (lower)
	r3 = mQ->p + x * mQ->step;
	for (int r = x; r < n; ++r, r3 += mQ->step) {
		memcpy(AUX2, r3 + x, (n - x) * sizeof(float));
		memset(r3, 0, n * sizeof(float));
		float *p3 = r3;
		for (int c = 0; c < n; ++c, ++p3) {
			float *p1 = AUX2, *p2 = Qgu + x * n + c;
			for (int k = x; k < n; ++k, ++p1, p2 += n)
				*p3 += *p1 * *p2;
		}
	}
}


static int deflate(const float *d, float *v, const int n, float *v1, float *d1, int *null_indexes, Rot_info *rotation_info) {
	// uses I2
	// return n1 (size of v1)
	int rot_numb = 0;
	// find identical entries in d and zero v elements with Givens rotation
	// We copy them temporarily to Z1
	memset(I2, 0, n * sizeof(int));  // elements to be skipped
	for (int i = 0; i < n; ++i) {
		if (I2[i]) continue;
		for (int j = i + 1; j < n; ++j) {
			if (fabsf(d[i] - d[j]) < threshold) {
				rotation_info[rot_numb].i = i;
				rotation_info[rot_numb].j = j;
				rotate(v[i], v[j], &rotation_info[rot_numb].c, &rotation_info[rot_numb].s);
				++rot_numb;
				v[i] = sqrtf(v[i]*v[i] + v[j]*v[j]);
				v[j] = 0;
				I2[j] = 1;
			}
		}
	}
	// terminator for rotation_info
	rotation_info[rot_numb].i = 0;
	rotation_info[rot_numb].j = 0;
	// Find Rotated v near 0 elements
	int n1 = 0, n2 = 0;
	for (int i = 0; i < n; ++i) {
		if (fabsf(v[i]) > threshold) {
			d1[n1] = d[i];
			v1[n1++] = v[i];
		}
		else null_indexes[n2++] = i;
	}
	null_indexes[n2] = -1;  // terminator
	for (int i = 0, j = n1; i < n2; ++i, ++j) d1[j] = d[null_indexes[i]];
	return n1;
}


static void gu_eisenstat(float *d, float *v, const int n, float *Qgu, float *dlam, const int Qgu_step) {
	// inputs: d, v (will be changed)
	// uses: I2, V
	const float norm_v = vector_norm(v, n);
	sort_indexes(I2, d, n);
	// sort v in Qgu (temporarily) and scale
	for (int i = 0; i < n; ++i) Qgu[i] = v[I2[i]] / norm_v;
	// sort d in v and scale
	for (int i = 0; i < n; ++i) v[i] = d[I2[i]] / (norm_v * norm_v);
	// copy v back to d
	memcpy(d, Qgu, n * sizeof(float));
	// now d and v are exchanged, so we swap the pointers
	float *tmp = d;
	d = v;
	v = tmp;
	// V = dvec (n*n)
	for (int k = 0; k < n; ++k) dlam[k] = zerodandc(d, v, n, k, V + k);
	for (int k = 0; k < n; ++k) {
		const float d_k_plus_1 = k < (n - 1) ? d[k + 1] : d[n - 1] + vector_norm_square(v, n);
		if (dlam[k] > 0) dlam[k] += d[k];
		else dlam[k] += d_k_plus_1;
	}
	for (int i = 0; i < n; ++i) {
		tmp = Qgu + i * Qgu_step;
		for (int j = 0; j < n; ++j) *tmp++ = 1;
	}
	for (int k = 0; k < n; ++k) {
		for (int old_j = 0; old_j < n; ++old_j) {
			int new_j = I2[old_j];
			if (fabsf(V[old_j * n + k]) > threshold) Qgu[new_j * Qgu_step + k] *= v[old_j] / V[old_j * n + k];
			else {
				float s = 0;
				for (int ii = 0; ii < n; ++ii) if (ii != old_j) s += v[ii] * v[ii] / V[ii * n + k];
				Qgu[new_j * Qgu_step + k] = -1 - s;
				for (int ii = 0; ii < n; ++ii) if (ii != new_j) Qgu[ii * Qgu_step + k] *= v[old_j];
			}
		}
		float column_norm = 0;
		tmp = Qgu + k;
		for (int ii = 0; ii < n; ++ii, tmp += Qgu_step) column_norm += *tmp * *tmp;
		column_norm = sqrtf(column_norm);
		tmp = Qgu + k;
		for (int ii = 0; ii < n; ++ii, tmp += Qgu_step) *tmp /= column_norm;
	}
	for (int k = 0; k < n; ++k) dlam[k] *= norm_v * norm_v;
}


#define SUM_VSQ_DL(A,B,R) do { R = 0; for (int ii = (A); ii < (B); ++ii) R += v[ii] * v[ii] / (d[ii] - lambda); } while(0)
#define SUM_VSQ_DLSQ(A,B,R) do { R = 0; for (int ii = (A); ii < (B); ++ii) R += v[ii] * v[ii] / ((d[ii] - lambda) * (d[ii] - lambda)); } while(0)


static float zerodandc(/*const*/ float *d, const float *v, const int n, const int i, float *dl) {
	// fills a column starting at dl
	// d is modified inside the function but restored before exiting
	float di = d[i];
	float di1, lambda;
	if (i < n - 1) {
		di1 = d[i + 1];
		lambda = (di + di1) / 2;
	}
	else {
		di1 = d[n - 1] +  vector_norm_square(v, n);
		lambda = di1;
	}
	float eta = 1, psi1, psi2;
	SUM_VSQ_DL(0, i + 1, psi1);
	SUM_VSQ_DL(i + 1, n, psi2);
	// NOTE: v is squared in original code, instead we square it when needed
	float offset;  // to restore d at the end
	if (1 + psi1 + psi2 > 0) {
		for (int ii = 0; ii < n; ++ii) d[ii] -= di;
		offset = di;
		lambda -= di;
		di1 -= di;
		while (fabsf(eta) > threshold) {
			float phi1, phi1s, psi2s;
			SUM_VSQ_DL(0, i, phi1);
			SUM_VSQ_DLSQ(0, i, phi1s);
			SUM_VSQ_DL(i + 1, n, psi2);
			SUM_VSQ_DLSQ(i + 1, n, psi2s);
			const float Di = -lambda, Di1 = di1 - lambda;
			const float a = (Di + Di1) * (1 + phi1 + psi2) - Di * Di1 * (phi1s + psi2s) + v[i] * v[i];
			const float b = Di * Di1 * (1 + phi1 + psi2) + Di1 * v[i] * v[i];
			const float c = (1 + phi1 + psi2) - Di * phi1s - Di1 * psi2s;
			if (a > 0) eta = (2 * b) /(a + sqrtf(a * a - 4 * b * c));
			else eta = (a - sqrtf(a * a - 4 * b * c)) / (2 * c);
			lambda += eta;
		}
	}
	else {
		for (int ii = 0; ii < n; ++ii) d[ii] -= di1;
		offset = di1;
		lambda -= di1;
		di -= di1;
		while (fabsf(eta) > threshold) {
			float psi1s, phi2, phi2s;
			SUM_VSQ_DL(0, i + 1, psi1);
			SUM_VSQ_DLSQ(0, i + 1, psi1s);
			SUM_VSQ_DL(i + 2, n, phi2);
			SUM_VSQ_DLSQ(i + 2, n, phi2s);
			const float Di = di - lambda, Di1 = - lambda;
			const float a = (Di + Di1) * (1 + psi1 + phi2) - Di * Di1 * (psi1s + phi2s) + v[i + 1] * v[i + 1];
			const float b = Di * Di1 * (1 + psi1 + phi2) + Di * v[i + 1] * v[i + 1];
			const float c = (1 + psi1 + phi2) - Di * psi1s - Di1 * phi2s;
			if (a > 0) eta = (2 * b) / (a + sqrtf(a * a - 4 * b * c));
			else eta = (a - sqrtf(a * a - 4 * b * c)) / (2 * c);
			lambda += eta;
		}
	}
	// fill column dl
	for (int ii = 0; ii < n; ++ii, dl += n) *dl = d[ii] - lambda;
	// restore d
	for (int ii = 0; ii < n; ++ii) d[ii] += offset;
	return lambda;
}


static const float *compare_values;
static int compare_indexes(const void *va, const void *vb) {
	const int *a = (const int*)va;
	const int *b = (const int*)vb;
	return (compare_values[*a] > compare_values[*b]) - (compare_values[*a] < compare_values[*b]);
}


static void sort_indexes(int *indexes, const float *values, int n) {
	for (int i = 0; i < n; ++i) indexes[i] = i;
	compare_values = values;
	qsort(indexes, n, sizeof(int), &compare_indexes);
}


