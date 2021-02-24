#include "alg.h"
#include <stdint.h>
#include <math.h>
#include <stm32f4xx.h>


volatile uint32_t tick;  // currently incremented every 1 ms

void timer_start(void) {
	tick = 0;
}

float timer_stop(void) {  // return ms
	static float period = 0;
	if (!period) period = (float)((SysTick->LOAD & 0xFFFFFF) + 1) / ((float)SystemCoreClock / 1000);  // ms
	return tick * period;
}

// precision time measurement using SysTick
// (measure very short times, limited by duration of single tick)
// Systick counts backwards from Reload to 0
// Counter and Reload are limited to 24 bits.
// Clock in this project is 180 MHz.
// This means: for Reload = 0xFFFFFF, max period is ~93.2 ms
//  for Reload = 180000 (as in this project), max duration is 1 ms
//  for Reload = 1800000, max period is 10 ms

static uint32_t precTick;

void systick_timer_start(void) {
	precTick = SysTick->VAL & 0xFFFFFF;
}

uint32_t systick_timer_stop(void) {
	return (precTick - (SysTick->VAL & 0xFFFFFF) + ((SysTick->LOAD & 0xFFFFFF) + 1)) % ((SysTick->LOAD & 0xFFFFFF) + 1);
}


float population[1];

float population_mean(unsigned count) {
	float sum = 0;
	for (float *p = population; p < (population + count); sum += *p++);
	return sum / count;
}

float population_stddev(unsigned count, float max_deviation_notused) {
	(void)max_deviation_notused;
	unsigned good = 0;
	float sum = 0, mean = population_mean(count);
	for (float *p = population; p < (population + count); ++p)
		if (1 /*fabsf(*p - mean) < max_deviation*/) {  // discard samples too far from mean
			sum += (*p - mean) * (*p - mean);
			++good;
		}
	return sqrtf(sum / good);
}


static float single_error(float n1, float n2) {
	// relative error
	return fabsf((n1 - n2) / n1);
}


float mean_error(float const *series1a, float const *series1b, unsigned count1, float const *series2a, float const *series2b, unsigned count2) {
	float e = 0;
	for (unsigned i = 0; i < count1; ++i) e += single_error(series1a[i], series1b[i]);
	if (series2a && series2b)
		for (unsigned i = 0; i < count2; ++i) e += single_error(series2a[i], series2b[i]);
	else count2 = 0;
	e /= count1 + count2;
	return e;
}


int compare_float_desc(const void *va, const void *vb) {
	// inspired from example in GNU libc documentation
	const float *a = (const float*)va;
	const float *b = (const float*)vb;
	return (*a < *b) - (*a > *b);
}
