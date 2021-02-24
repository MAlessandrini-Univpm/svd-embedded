#include "c_math.h"
#include <string.h>
#include <math.h>
#include <assert.h>

extern float AUX1[], P[], Q[];

//#define STUPID_MUL  // don't use optimized multiplication for special matrices, for comparisons in the paper


void sub_matrix (float *d, const float *s, int rows, int cols, int src_step, bool transpose) {
	// copy (and possibly transpose) top-left portion of s into d
	// d and s must not overlap!
	// s is row-major ordering!
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			if (transpose)
				d[j * rows + i] = s[i * src_step + j];
			else
				d[i * cols + j] = s[i * src_step + j];
}


float sign(float n) {
	if (n == 0) return 0;
	if (n > 0) return 1;
	return -1;
}


float signNo0(float n) {
	if (n < 0) return -1;
	return 1;
}


void eye(float *m, int n) {
	memset(m, 0, n * n * sizeof(float));
	for (int i = 0; i < n; ++i, m += n + 1) *m = 1;
}


void householder(float *m, int rows, int columns, float *U, float *V) {
	// uses AUX1, P, Q
	static const float thresh = 1e-7;
	if (U) eye(U, rows);
	if (V) eye(V, columns);
	for (int k = 0; k < columns; ++k) {
		float *v = AUX1;
		float *p = m + k * columns + k;
		for (int i = k; i < rows; ++i, p += columns) *v++ = *p;
		v = AUX1;
		int sv = rows - k;  // size of v
		assert(sv > 0);
		//for (int i = 0; i < sv; ++i) assert(!isnan(v[i]) && !isinf(v[i]));
		float norm = vector_norm(v, sv);
		//assert(!isnan(norm) && !isinf(norm));
		v[0] += signNo0(v[0]) * norm;
		norm = vector_norm(v, sv);
		//assert(!isnan(norm) && !isinf(norm));
		if (norm > thresh) for (int i = 0; i < sv; ++i) v[i] /= norm;
		// Q will hold 2 * v * v' [sv x sv]
		p = Q;
		for (int i = 0; i < sv; ++i)
			for (int j = 0; j < sv; ++j)
				*p++ = 2 * v[i] * v[j];
		// P will hold the product
		Mat m1 = { Q, sv, sv, sv }, m2 = { m + k * columns + k, sv, columns - k, columns },
			m3 = { P, sv, columns - k, columns - k };
		matrix_mul(&m1, &m2, &m3);
		for (int i = k; i < rows; ++i) {
			p = m + i * columns + k;
			for (int j = k; j < columns; ++j)
				*p++ -= *m3.p++;
		}
		if (U) {  // U is rows*rows
			Mat m1 = { U + k, rows, sv, rows }, m2 = { Q, sv, sv, sv },
				m3 = { P, rows, sv, sv };
			matrix_mul(&m1, &m2, &m3);
			for (int i = 0; i < rows; ++i) {
				p = U + i * rows + k;
				for (int j = k; j < rows; ++j)
					*p++ -= *m3.p++;
			}
		}
		if (k <= columns - 3) {
			sv = columns - k - 1;
			v = AUX1;
			memcpy(v, m + k * columns + k + 1, sv * sizeof(float));
			norm = vector_norm(v, sv);
			v[0] += signNo0(v[0]) * norm;
			norm = vector_norm(v, sv);
			if (norm > thresh) for (int i = 0; i < sv; ++i) v[i] /= norm;
			p = Q;
			for (int i = 0; i < sv; ++i)
				for (int j = 0; j < sv; ++j)
					*p++ = 2 * v[i] * v[j];
			Mat m1 = { m + k * columns + k + 1, rows - k, sv, columns }, m2 = { Q, sv, sv, sv },
				m3 = { P, rows - k, sv, sv };
			matrix_mul(&m1, &m2, &m3);
			for (int i = k; i < rows; ++i) {
				p = m + i * columns + k + 1;
				for (int j = k + 1; j < columns; ++j)
					*p++ -= *m3.p++;
			}
			if (V) {  // v is columns*columns
				Mat m1 = { V + k + 1, columns, sv, columns }, m2 = { Q, sv, sv, sv },
					m3 = { P, columns, sv, sv };
				matrix_mul(&m1, &m2, &m3);
				for (int i = 0; i < columns; ++i) {
					p = V + i * columns + k + 1;
					for (int j = k + 1; j < columns; ++j)
						*p++ -= *m3.p++;
				}
			}
		}
	}
}


float vector_norm(const float *v, int n) {
	float norm = 0;
	while (n--) {
		//assert(!isnan(*v) && !isinf(*v));
		//assert(fabsf(*v) < 1e20);
		norm += *v * *v;
		//assert(!isnan(norm) && !isinf(norm));
		++v;
	}
	return sqrtf(norm);
}


float vector_norm_square(const float *v, int n) {
	float norm = 0;
	while (n--) {
		norm += *v * *v;
		++v;
	}
	return norm;
}


void matrix_mul(const Mat *m1, const Mat *m2, Mat *m3) {
	// m3 = m1 * m2
	const float *r1 = m1->p;
	float *r3 = m3->p;
	for (int r = 0; r < m3->rows; ++r, r3 += m3->step, r1 += m1->step) {
		memset(r3, 0, m3->columns * sizeof(float));
		float *p3 = r3;
		for (int c = 0; c < m3->columns; ++c, ++p3) {
			const float *p1 = r1, *p2 = m2->p + c;
			for (int k = 0; k < m1->columns; ++k, ++p1, p2 += m2->step)
				*p3 += *p1 * *p2;
		}
	}
}


void matrix_mul_symm(const Mat *m1, Mat *m3) {
	// m3 = m1 * m1'
	// result is symmetric matrix
	// matrices must be square
#ifdef STUPID_MUL
	// full matrix mul, with m2 replaced by m1 swept in transposed way
	const float *r1 = m1->p;
	float *r3 = m3->p;
	for (int r = 0; r < m3->rows; ++r, r3 += m3->step, r1 += m1->step) {
		memset(r3, 0, m3->columns * sizeof(float));
		float *p3 = r3;
		for (int c = 0; c < m3->columns; ++c, ++p3) {
			const float *p1 = r1, *p2 = m1->p + c * m1->step;
			for (int k = 0; k < m1->columns; ++k, ++p1, ++p2)
				*p3 += *p1 * *p2;
		}
	}
#else
	const float *r1 = m1->p;
	float *r3 = m3->p;
	// compute only half
	for (int r = 0; r < m3->rows; ++r, r3 += m3->step, r1 += m1->step) {
		memset(r3, 0, (r + 1) * sizeof(float));
		float *p3 = r3;
		for (int c = 0; c <= r ; ++c, ++p3) {  // only half columns
			const float *p1 = r1, *p2 = m1->p + c * m1->step;
			for (int k = 0; k < m1->rows; ++k)
				*p3 += *p1++ * *p2++;
		}
	}
	// fill other half
	r3 = m3->p + 1;  // start of first half-line
	const int diff = m3->step - m3->columns;
	for (int r = 0; r < m3->rows - 1; ++r) {
		r1 = r3 + m3->step - 1;
		for (int c = 0; c < m3->rows - r - 1; ++c, r1 += m3->step) *r3++ = *r1;
		r3 += diff + r + 2;  // set r3 to start of next half-line
	}
#endif
}


void matrix_mul_symm_bidiag(const Mat *m1, Mat *m3) {
	// m3 = m1 * m1'  with m1 upper bidiagonal
	// result is symmetric matrix
	// matrices must be square
#ifdef STUPID_MUL
	matrix_mul_symm(m1, m3);
#else
	const float *r1 = m1->p;
	float *r3 = m3->p;
	// result is tridiagonal, so most elements will be zero
	for (int k = 0; k < m3->rows; ++k, r3 += m3->step) memset(r3, 0, m3->columns * sizeof(float));
	r3 = m3->p;
	// compute only half
	// upper-left element:
	*r3 = *r1 * *r1;
	++r1;
	*r3 += *r1 * *r1;
	r1 += m1->step;
	r3 += m3->step;
	// rows except first and last one
	for (int count = 0; count < m3->rows - 2; ++count, r3 += m3->step, r1 += m1->step) {
		*r3++ = *r1 * *(r1 - m1->step);
		*r3 = *r1 * *r1;
		++r1;
		*r3 += *r1 * *r1;
	}
	// last row
	*r3++ = *r1 * *(r1 - m1->step);
	*r3 = *r1 * *r1;
	// fill other half
	r3 = m3->p + 1;  // first element of upper diagonal
	for (int r = 0; r < m3->rows - 1; ++r, r3 += m3->step + 1)
		*r3 = *(r3 + m3->step - 1);  // matching element of lower diagonal
#endif
}


void givens_rotation_left(float c, float s, int i, int k, bool tran, Mat *mat) {
	// perform mat = G * mat, multiplying only the two rows affected
	const float e = tran ? -1 : 1;
	float *r1 = mat->p + i * mat->step, *r2 = mat->p + k * mat->step;
	for (int kk = 0; kk < mat->columns; ++kk, ++r1, ++r2) {
		const float temp = *r1;
		*r1 = c * *r1 + s * e * *r2;
		*r2 = -s * e * temp + c * *r2;
	}
}


void givens_rotation_right(float c, float s, int i, int k, bool tran, Mat *mat) {
	// perform mat = mat * G, multiplying only the two columns affected
	const float e = tran ? -1 : 1;
	float *c1 = mat->p + i, *c2 = mat->p + k;
	for (int kk = 0; kk < mat->rows; ++kk, c1 += mat->step, c2 += mat->step) {
		const float temp = *c1;
		*c1 = c * *c1 - s * e * *c2;
		*c2 = s * e * temp + c * *c2;
	}
}


void golub_kahan_step(float *B, const int n, const int p, const int q, float *U, float *V, const int rows) {
	// uses P
	const int dim = n - q - p - 2;
	// P = product
	Mat m1 = { B + (p + 1) * n + p + 1, dim, dim, n }, mC = { P, dim, dim, dim };
	matrix_mul_symm_bidiag(&m1, &mC);
	mC.p += (dim - 2) * dim + dim - 2;
	mC.rows = 2;
	mC.columns = 2;
	float lambda_1,lambda_2;
	eig2(&mC, &lambda_1, &lambda_2);
	const float mu = fabsf(mC.p[mC.step + 1] - lambda_1) < fabsf(mC.p[mC.step + 1] - lambda_2) ? lambda_1 : lambda_2;
	const int k = p + 1;
	float alfa = B[k * n + k] * B[k * n + k] - mu;
	float beta  = B[k * n + k] * B[k * n + k + 1];
	for (int k = p + 1; k <= n - q - 3; ++k) {
		float c, s;
		rotate(alfa, beta, &c, &s);
		Mat m = { B, n, n, n };
		givens_rotation_right(c, s, k, k + 1, true, &m);
		if (V) {
			m.p = V;
			givens_rotation_right(c, s, k, k + 1, true, &m);
		}
		alfa = B[k * n + k];
		beta  = B[(k + 1) * n + k];
		rotate(alfa, beta, &c, &s);
		m.p = B;
		givens_rotation_left(c, s, k, k + 1, false, &m);
		if (U) {
			m.p = U;
			m.rows = rows;
			m.step = rows;
			givens_rotation_right(c, s, k, k + 1, true, &m);
		}
		alfa = B[k * n + k + 1];
		if (k <= (n - q - 4)) beta = B[k * n + k + 2];
	}
}


void eig2(const Mat *m, float *l1, float *l2) {
	const float b = m->p[m->step + 1] + m->p[0];
	const float r = sqrtf(b * b - 4 * (m->p[0] * m->p[m->step + 1] - m->p[m->step] * m->p[1]));
	*l1 = (b + r) / 2;
	*l2 = (b - r) / 2;
}


void rotate(float a, float b, float *c, float *s) {
	static const float thresh = 1e-7;
	if (fabsf(b) < thresh) {
		*c = 1;
		*s = 0;
	}
	else if (fabsf(a) < thresh) {
		*c = 0;
		*s = sign(b);
	}
	else {
		*c = a / sqrtf(a * a + b * b);
		*s = b / sqrtf(a * a + b * b);
	}
}
