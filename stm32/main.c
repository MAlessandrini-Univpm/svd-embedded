#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "alg.h"
#include "c_math.h"
#include "stm32f4xx_conf.h"


//#define RATIO_2  // matrix ratio of 1:2 instead of 3:4
//#define RATIO_1  // matrix ratio of 1:1

// max matrix dimension
// NOTE: MAXDIM can be up to 144 removing U, V and the call to svd_divide_conquer_C
#ifdef RATIO_2
#define DIM_DIV 2
#define DIM_MUL 1
#elif defined RATIO_1
#define DIM_DIV 1
#define DIM_MUL 1
#else
#define DIM_DIV 4
#define DIM_MUL 3
#endif

#define MAXDIM2 72  // max columns, keep those constant
#define MAXDIM (MAXDIM2 / DIM_MUL * DIM_DIV)  // max rows, depending on ratio


// define this macro if measuring energy consumption
// with an external device -> no debug, pauses between measures, etc.
//#define EXTERNAL_MEASURE


static const float A[200 * 150] = {  // should be placed in ROM
#include "A200x150.h"
};

float B[MAXDIM * MAXDIM2] __attribute__((section(".ccmram")));
float P[MAXDIM * MAXDIM], Q[MAXDIM * MAXDIM], Z1[MAXDIM], AUX1[MAXDIM] __attribute__((section(".ccmram"))), AUX2[MAXDIM] __attribute__((section(".ccmram"))), Z3[MAXDIM] __attribute__((section(".ccmram")));
float U[MAXDIM2 * MAXDIM2] __attribute__((section(".ccmram"))), V[MAXDIM2 * MAXDIM2] , UU[0];  // U,V only used by D.C.
int I1[MAXDIM], I2[MAXDIM];
Rot_info Rot[MAXDIM];
uint32_t precisionTimer = 0;  // number of SysTick increments

// matlab reference SVD results, should be placed in ROM
// ratio 3:4
static const float B32x24[] = {
#include "ref_matlab/B32x24"
};
static const float B48x36[] = {
#include "ref_matlab/B48x36"
};
static const float B64x48[] = {
#include "ref_matlab/B64x48"
};
static const float B80x60[] = {
#include "ref_matlab/B80x60"
};
static const float B96x72[] = {
#include "ref_matlab/B96x72"
};
// ratio 1:2
static const float B48x24[] = {
#include "ref_matlab/B48x24"
};
static const float B72x36[] = {
#include "ref_matlab/B72x36"
};
static const float B96x48[] = {
#include "ref_matlab/B96x48"
};
static const float B120x60[] = {
#include "ref_matlab/B120x60"
};
static const float B144x72[] = {
#include "ref_matlab/B144x72"
};
// ratio 1:1
static const float B24x24[] = {
#include "ref_matlab/B24x24"
};
static const float B36x36[] = {
#include "ref_matlab/B36x36"
};
static const float B48x48[] = {
#include "ref_matlab/B48x48"
};
static const float B60x60[] = {
#include "ref_matlab/B60x60"
};
static const float B72x72[] = {
#include "ref_matlab/B72x72"
};


void __attribute__ ((noinline)) exit_with_error(volatile int line) {
	(void)line;  // to be inspected from debugger
	while(1) {}
}

int __errno;  // needed by sqrtf
void _exit() {  // needed somewhere in C runtime (not with --specs=nosys.specs, but this is more useful)
	exit_with_error(__LINE__);
}


static void do_pause() {
#ifdef EXTERNAL_MEASURE
	extern volatile uint32_t tick;
	tick = 0;
	while (tick < 1000) __WFI();
#endif
}

void setGpios();

const int rows_A = 200, cols_A = 150;

static int single_test(int rows, int cols, int(*f)(int,int,float*,float*), float *mean, float *stddev, unsigned iterations_notused, float max_deviation, bool transpose) {
	// return iterations performed by algorithm (not to be confused with number of calls that was in Raspberry code)
	(void)iterations_notused;
	sub_matrix(B, A, rows, cols, cols_A, transpose); // some algorithms want column major
	timer_start();
	int alg_iter = f(rows, cols, NULL, NULL);  // call actual function
	population[0] = timer_stop();
	*mean = population_mean(1);
	*stddev = population_stddev(1, max_deviation);
	return alg_iter;
}

// we declare table_line volatile to prevent the compiler to optimize it out,
// so we can inspect it in the debugger
static volatile float table_line[15];  // to output a CSV line

static const int colsSweep[] = { 24, 36, 48, 60, 72, 0 };  // null-terminated


int main () {
	if (SysTick_Config(SystemCoreClock / 1000)) exit_with_error(__LINE__);  // configure systick for interrupt every 1 ms (increment time variable every 1 ms)

#ifdef EXTERNAL_MEASURE
	setGpios();
#endif

	while(1) {
		for (const int *sweep = colsSweep; *sweep; ++sweep) {
			const int cols_sub = *sweep, rows_sub = cols_sub / DIM_MUL * DIM_DIV;

			// reference vector
			const float *Z2;
			if (cols_sub == 24 && rows_sub == 32) Z2 = B32x24;
			else if (cols_sub == 36 && rows_sub == 48) Z2 = B48x36;
			else if (cols_sub == 48 && rows_sub == 64) Z2 = B64x48;
			else if (cols_sub == 60 && rows_sub == 80) Z2 = B80x60;
			else if (cols_sub == 72 && rows_sub == 96) Z2 = B96x72;
			else if (cols_sub == 24 && rows_sub == 48) Z2 = B48x24;
			else if (cols_sub == 36 && rows_sub == 72) Z2 = B72x36;
			else if (cols_sub == 48 && rows_sub == 96) Z2 = B96x48;
			else if (cols_sub == 60 && rows_sub == 120) Z2 = B120x60;
			else if (cols_sub == 72 && rows_sub == 144) Z2 = B144x72;
			else if (cols_sub == 24 && rows_sub == 24) Z2 = B24x24;
			else if (cols_sub == 36 && rows_sub == 36) Z2 = B36x36;
			else if (cols_sub == 48 && rows_sub == 48) Z2 = B48x48;
			else if (cols_sub == 60 && rows_sub == 60) Z2 = B60x60;
			else if (cols_sub == 72 && rows_sub == 72) Z2 = B72x72;
			else exit_with_error(__LINE__);

			table_line[0] = cols_sub;
			table_line[1] = rows_sub;  // changed from previous code

			float mean1 = 0, stddev1 = 0, e;

#ifdef EXTERNAL_MEASURE
			const int repeat = cols_sub <= 24 ? 100 : cols_sub <= 36 ? 40 : cols_sub <= 48 ? 20 : 5;
#else
			const int repeat = 1;
#endif

			do_pause();

			// Golub-Reinsch, implementation in C
			for (int rep = 0; rep < repeat; ++rep) single_test(rows_sub, cols_sub, &svd_golub_reinsch_C, &mean1, &stddev1, 1, 0, false);
			for (int i = 0; i < cols_sub; ++i) AUX1[i] = B[i * cols_sub + i];
			qsort(AUX1, cols_sub, sizeof(float), &compare_float_desc);
			e = mean_error(Z2, AUX1, cols_sub, NULL, NULL, 0);
			table_line[5] = mean1;
			table_line[6] = e;

			do_pause();

			// 1-sided Jacobi, implementation in C
			for (int rep = 0; rep < repeat; ++rep) single_test(rows_sub, cols_sub, &svd_one_sided_jacobi_C, &mean1, &stddev1, 1, 0, !ONE_SIDED_JACOBI_ROW_MAJOR);
			qsort(AUX1, cols_sub, sizeof(float), &compare_float_desc);
			e = mean_error(Z2, AUX1, cols_sub, NULL, NULL, 0);
			table_line[7] = mean1;
			table_line[8] = e;

			do_pause();

			// Jacobi rotation, implementation in C
			for (int rep = 0; rep < repeat; ++rep) single_test(rows_sub, cols_sub, &svd_jacobi_rotation_C, &mean1, &stddev1, 1, 0, false);
			qsort(AUX1, cols_sub, sizeof(float), &compare_float_desc);
			e = mean_error(Z2, AUX1, cols_sub, NULL, NULL, 0);
			table_line[9] = mean1;
			table_line[10] = e;

			do_pause();

			// Demmel-Kahan, implementation in C
			for (int rep = 0; rep < repeat; ++rep) single_test(rows_sub, cols_sub, &svd_demmel_kahan_C, &mean1, &stddev1, 1, 0, false);
			for (int i = 0; i < cols_sub; ++i) AUX1[i] = B[i * cols_sub + i];
			qsort(AUX1, cols_sub, sizeof(float), &compare_float_desc);
			e = mean_error(Z2, AUX1, cols_sub, NULL, NULL, 0);
			table_line[11] = mean1;
			table_line[12] = e;

			do_pause();

			// Divide and Conquer
			for (int rep = 0; rep < repeat; ++rep) single_test(rows_sub, cols_sub, &svd_divide_conquer_C, &mean1, &stddev1, 1, 0, false);
			qsort(AUX1, cols_sub, sizeof(float), &compare_float_desc);
			e = mean_error(Z2, AUX1, cols_sub, NULL, NULL, 0);
			table_line[13] = mean1;
			table_line[14] = e;

			//table_line[2] = (float)(precisionTimer) / ((float)SystemCoreClock / 1000);  // ms

#ifndef EXTERNAL_MEASURE
			__asm__ __volatile__ ("bkpt #0");
#else
			do_pause(); do_pause(); do_pause();
#endif
		}
#ifdef EXTERNAL_MEASURE
		while(1) __WFI();
#endif
	}
}




void setGpios() {
	// set pulldown to all GPIOs, except pins with special default settings,
	//  so to minimize power consumption
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOA, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOB, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOC, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOD, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOE, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOF, ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOG, ENABLE);
	GPIO_InitTypeDef gpio;
	GPIO_StructInit(&gpio);  // all pins, input, Output pull-push, 2 MHz, no PU/PD
	gpio.GPIO_PuPd = GPIO_PuPd_DOWN;
	GPIO_Init(GPIOC, &gpio);
	GPIO_Init(GPIOD, &gpio);
	GPIO_Init(GPIOE, &gpio);
	GPIO_Init(GPIOF, &gpio);
	GPIO_Init(GPIOG, &gpio);
	gpio.GPIO_Pin = GPIO_Pin_All & ~(GPIO_Pin_13 | GPIO_Pin_14 | GPIO_Pin_15);
	GPIO_Init(GPIOA, &gpio);
	gpio.GPIO_Pin = GPIO_Pin_All & ~(GPIO_Pin_3 | GPIO_Pin_4);
	GPIO_Init(GPIOB, &gpio);
}
