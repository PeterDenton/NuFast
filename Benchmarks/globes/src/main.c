#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "glb_probability.h"

#define M_PI 3.141592653589793238462643383279

// gets the difference in time in nanoseconds
int difftimespec_ns(const struct timespec end_time, const struct timespec start_time)
{
	return 1000 * (((int)end_time.tv_sec - (int)start_time.tv_sec) * (int)1000000 + ((int)end_time.tv_nsec - (int)start_time.tv_nsec) / 1000);
}
double Speed(int n)
{
	/* Setup default probability engine */
	glb_init_probability_engine();

	osc_params_angles opa;
	opa.theta12 = 33 * M_PI / 180;
	opa.theta13 = 8.5 * M_PI / 180;
	opa.theta23 = 48 * M_PI / 180;
	opa.delta = 0.7 * M_PI;
	opa.dms = 7.5e-5;
	opa.dma = 2.5e-3;
	glb_set_oscillation_parameters(opa);

	double P[3][3];
	int cp_sign = +1;
	double rho = 3.;
	double L = 1300.;

	// initialize the grid
	double *Es = malloc(n * sizeof(double));

	double E_min, E_scale, E_max;

	E_min = 1e-1;
	E_max = 1e2;

	E_scale = pow(E_max / E_min, 1. / n);

	// fill in arrays outside of timing
	for (int i = 0; i < n; i++)
		Es[i] = E_min * pow(E_scale, i);

	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
	for (int i = 0; i < n; i++)
		glb_probability_matrix(P, cp_sign, Es[i], &L, &rho);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);

	glb_free_probability_engine();
	free(Es);

	return 1.0 * difftimespec_ns(end_time, start_time) / n;
}

int main()
{
	double s, speed_sum, speedsq_sum, mean, std;
	int m, n;

	m = 1e3;
	n = 1e5;

	speed_sum = 0;
	speedsq_sum = 0;
	for (int i = 0; i < m; i++)
	{
if (i % 100 == 0) printf("%g\n", 1.0 * i / m);
		s = Speed(n);
		speed_sum += s;
		speedsq_sum += s * s;
	} // i, m
	mean = speed_sum / m;
	std = sqrt(speedsq_sum / m - mean * mean);
	printf("t = %g +- %g ns\n", mean, std);

	return 0;
}
