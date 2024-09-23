#include <cmath>
#include <cstdio>
#include <chrono>
#include "Benchmark.h"
#include "NuFast.h"
#include "Exact_Cubic.h"
#include "Page.h"
#include "DMP.h"

// Some constants
double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
double const YerhoE2a = 1.52588e-4;
double const oneThird = 1. / 3;

// Compute the speed n times and calculate the mean over those n times
// N_Newton is used for NuFast only, ignored otherwise
double Speed_Helper(void Probability_Calculator(double, double, double, double, double, double, double, double, double, double, int, double (*probs_returned)[3][3]), int n, int N_Newton)
{
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, Emin, Emax, Estep, rho, Ye, probs_returned[3][3];

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = -0.7 * M_PI;
	Dmsq21 = 7.5e-5; // eV^2
	Dmsq31 = 2.5e-3; // eV^2

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	Emin = 0.5; // GeV
	Emax = 5;
	Estep = (Emax - Emin) / n;
	rho = 3; // g/cc
	Ye = 0.5;

	std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n; i++)
	{
		E = Emin + i * Estep;
		Probability_Calculator(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);
	} // i, E, n
	std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / n;
}

// Compute the speed many times and calculate the mean and standard deviation
// m is outer loop, n, is inner loop
// N_Newton is used for NuFast only, ignored otherwise
void Speed(void Probability_Calculator(double, double, double, double, double, double, double, double, double, double, int, double (*probs_returned)[3][3]), int m, int n, int N_Newton, double &mean, double &std)
{
	double s, speed_sum, speedsq_sum;

	speed_sum = 0;
	speedsq_sum = 0;
	for (int i = 0; i < m; i++)
	{
		s = Speed_Helper(Probability_Calculator, n, N_Newton) * 1e9; // in ns
		speed_sum += s;
		speedsq_sum += s * s;
	} // i, m
	mean = speed_sum / m;
	std = sqrt(speedsq_sum / m - mean * mean);
}

int main()
{
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
	double L, E, rho, Ye;
	double Emin, Emax, Estep;
	double mean, std;
	int N_Newton, m, n;
	double probs_returned[3][3], probs_returned_Exact_Cubic[3][3];

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = -0.7 * M_PI;
	Dmsq21 = 7.5e-5;  // eV^2
	Dmsq31 = +2.5e-3; // eV^2

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	E = 2.5; // GeV
	rho = 3; // g/cc
	Ye = 0.5;

	N_Newton = 0; // <-- Sets the precision. 0 is close to the single precision limit and is better than DUNE/HK in the high statistics regime. Increasig N_Newton improves the precision at a modest computational cost

	// ------------------------------------------ //
	// Calculate all 9 oscillationa probabilities //
	// ------------------------------------------ //
	Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);
	Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &probs_returned_Exact_Cubic);

	// ------------------------------- //
	// Print out the fractional errors //
	// ------------------------------- //
	printf("L = %g E = %g rho = %g Ye = %g N_Newton = %d\n", L, E, rho, Ye, N_Newton);
	printf("Precision (N_Newton = %d):\n", N_Newton);
	printf("alpha beta DeltaP/P\n");
	for (int alpha = 0; alpha < 3; alpha++)
	{
		for (int beta = 0; beta < 3; beta++)
		{
			printf("%d %d %g\n", alpha, beta, fabs(probs_returned[alpha][beta] - probs_returned_Exact_Cubic[alpha][beta]) / probs_returned_Exact_Cubic[alpha][beta]);
		} // beta, 3
	} // alpha, 3

	// ----------------------------------- //
	// Write the fractional errors to file //
	// ----------------------------------- //
	n = 1e3;
	// DUNE //
	Emin = 0.5;
	Emax = 5;
	Estep = (Emax - Emin) / n;
	FILE *dataf = fopen("data/Precision_DUNE.txt", "w");
	fprintf(dataf, "%g %g %g\n", L, rho, Ye);
	for (int i = 0; i <= n; i++)
	{
		E = Emin + i * Estep;
		fprintf(dataf, "%g ", E);
		Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &probs_returned_Exact_Cubic);
		for (int j = 0; j < 2; j++)
		{
			Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, j, &probs_returned);
			fprintf(dataf, "%g %g ", fabs(probs_returned[1][1] - probs_returned_Exact_Cubic[1][1]) / probs_returned_Exact_Cubic[1][1], fabs(probs_returned[1][0] - probs_returned_Exact_Cubic[1][0]) / probs_returned_Exact_Cubic[1][0]);
		} // j, N_Newton, 2
		// now do DMP
		Probability_Matter_LBL_DMP(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &probs_returned);
		fprintf(dataf, "%g %g ", fabs(probs_returned[1][1] - probs_returned_Exact_Cubic[1][1]) / probs_returned_Exact_Cubic[1][1], fabs(probs_returned[1][0] - probs_returned_Exact_Cubic[1][0]) / probs_returned_Exact_Cubic[1][0]);
		fprintf(dataf, "\n");
	} // i, E, n
	fclose(dataf);

	// HK //
	L = 295;
	Emin = 0.1;
	Emax = 2;
	Estep = (Emax - Emin) / n;
	dataf = fopen("data/Precision_HK.txt", "w");
	fprintf(dataf, "%g %g %g\n", L, rho, Ye);
	for (int i = 0; i <= n; i++)
	{
		E = Emin + i * Estep;
		fprintf(dataf, "%g ", E);
		Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &probs_returned_Exact_Cubic);
		for (int j = 0; j < 2; j++)
		{
			Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, j, &probs_returned);
			fprintf(dataf, "%g %g ", fabs(probs_returned[1][1] - probs_returned_Exact_Cubic[1][1]) / probs_returned_Exact_Cubic[1][1], fabs(probs_returned[1][0] - probs_returned_Exact_Cubic[1][0]) / probs_returned_Exact_Cubic[1][0]);
		} // j, N_Newton, 2
		// now do DMP
		Probability_Matter_LBL_DMP(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, 0, &probs_returned);
		fprintf(dataf, "%g %g ", fabs(probs_returned[1][1] - probs_returned_Exact_Cubic[1][1]) / probs_returned_Exact_Cubic[1][1], fabs(probs_returned[1][0] - probs_returned_Exact_Cubic[1][0]) / probs_returned_Exact_Cubic[1][0]);
		fprintf(dataf, "\n");
	} // i, E, n
	fclose(dataf);

	// JUNO //
	L = 50;
	rho = 2.6;
	Emin = 1e-3;
	Emax = 10e-3;
	Estep = (Emax - Emin) / n;
	dataf = fopen("data/Precision_JUNO.txt", "w");
	fprintf(dataf, "%g %g %g %d\n", L, rho, Ye, N_Newton);
	for (int i = 0; i <= n; i++)
	{
		E = Emin + i * Estep;
		Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, -E, rho, Ye, N_Newton, &probs_returned);
		Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, -E, rho, Ye, 0, &probs_returned_Exact_Cubic);
		fprintf(dataf, "%g %g\n", E, fabs(probs_returned[0][0] - probs_returned_Exact_Cubic[0][0]) / probs_returned_Exact_Cubic[0][0]);
	} // i, E, n
	fclose(dataf);

	// ---------- //
	// Speed test //
	// ---------- //
	dataf = fopen("data/Speed.txt", "w");
	m = int(1e3);
	n = int(1e6);

	printf("\nSpeed (this takes a moment):\n");

	// Vacuum probability
	printf("Vacuum (1/7)\n");
	Speed(Probability_Vacuum_LBL, m, n, 0, mean, std);
	printf("%g +- %g ns\n", mean, std);
	fprintf(dataf, "Vacuum %g += %g ns\n", mean, std);

	// NuFast at various N_Newton levels
	for (N_Newton = 0; N_Newton <= 3; N_Newton++)
	{
		printf("N_Newton = %d (%d/7)\n", N_Newton, N_Newton + 2);
		Speed(Probability_Matter_LBL, m, n, N_Newton, mean, std);
		printf("%g +- %g ns\n", mean, std);
		fprintf(dataf, "N_Newton=%d %g += %g ns\n", N_Newton, mean, std);
	} // N_Newton, 0, 3

	// Exact with cubic from Cardano/ZS
	printf("Exact_Cubic (6/7)\n");
	Speed(Probability_Matter_LBL_Exact_Cubic, m, n, 0, mean, std);
	printf("%g +- %g ns\n", mean, std);
	fprintf(dataf, "Exact_Cubic %g += %g ns\n", mean, std);

	// Page's algorithm from 2309.06900
	printf("Page (7/7)\n");
	Speed(Probability_Matter_LBL_Page, m, n, 0, mean, std);
	printf("%g +- %g ns\n", mean, std);
	fprintf(dataf, "Page %g += %g ns\n", mean, std);

	fclose(dataf);

	return 0;
}
