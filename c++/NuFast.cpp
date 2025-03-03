#include <cmath>
#include <stdio.h>

// Some constants
constexpr double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
constexpr double const YerhoE2a = 1.52588e-4;

// Probability_Matter_LBL calculates all nine oscillation probabilities including
// the matter effect in an optimized, fast, and efficient way. The precision can
// be controlled with N_Newton. For many applications N_Newton=0 may be enough,
// but many years of DUNE or HK-LBL may require N_Newton=1. This code may be
// suitable for atmospheric neutrinos. The code is standalone.
//
// Inputs:
//   mixing angles (usual parameterization)
//   phase (usual parameterization) make Dmsq31 positive/negative for the NO/IO
//   Delta msq's (eV^2)
//   L (km)
//   E (GeV) positive for neutrinos, negative for antineutrinos
//   rho (g/cc)
//   Ye: electron fraction, typically around 0.5
//   N_Newton: number of Newton's method iterations to do. should be zero, one, two (or higher)
// Outputs:
//   probs_returned is all nine oscillation probabilities: e.g. probs_returned[1][0] is mu->e
void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, double E, double rho, double Ye, int N_Newton, double (*probs_returned)[3][3])
{
	double c13sq, sind, cosd, Jrr, Jmatter, Dmsqee, Amatter;
	double Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
	double A, B, C;
	double See, Tee, Smm, Tmm;
	double xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
	double Xp2, Xp3, PiDlambdaInv, tmp;
	double Lover4E, D21, D32;
	double sinD21, sinD31, sinD32;
	double sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
	double Pme_CPC, Pme_CPV, Pmm, Pee;

	// --------------------------------------------------------------------- //
	// First calculate useful simple functions of the oscillation parameters //
	// --------------------------------------------------------------------- //
	c13sq = 1 - s13sq;

	// Ueisq's
	Ue2sq = c13sq * s12sq;
	Ue3sq = s13sq;

	// Umisq's, Utisq's and Jvac	 
	Um3sq = c13sq * s23sq;
	// Um2sq and Ut2sq are used here as temporary variables, will be properly defined later	 
	Ut2sq = s13sq * s12sq * s23sq;
	Um2sq = (1 - s12sq) * (1 - s23sq);

	Jrr = sqrt(Um2sq * Ut2sq);
	sind = sin(delta);
	cosd = cos(delta);

	Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd;
	Jmatter = 8 * Jrr * c13sq * sind;
	Amatter = Ye * rho * E * YerhoE2a;
	Dmsqee = Dmsq31 - s12sq * Dmsq21;

	// calculate A, B, C, See, Tee, and part of Tmm
	A = Dmsq21 + Dmsq31; // temporary variable
	See = A - Dmsq21 * Ue2sq - Dmsq31 * Ue3sq;
	Tmm = Dmsq21 * Dmsq31; // using Tmm as a temporary variable	  
	Tee = Tmm * (1 - Ue3sq - Ue2sq);
	C = Amatter * Tee;
	A = A + Amatter;

	// ---------------------------------- //
	// Get lambda3 from lambda+ of MP/DMP //
	// ---------------------------------- //
	xmat = Amatter / Dmsqee;
	tmp = 1 - xmat;
	lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + sqrt(tmp * tmp + 4 * s13sq * xmat));

	// ---------------------------------------------------------------------------- //
	// Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) //
	// ---------------------------------------------------------------------------- //
	B = Tmm + Amatter * See; // B is only needed for N_Newton >= 1
	for (int i = 0; i < N_Newton; i++)
		lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B); // this strange form prefers additions to multiplications

	// ------------------- //
	// Get  Delta lambda's //
	// ------------------- //
	tmp = A - lambda3;
	Dlambda21 = sqrt(tmp * tmp - 4 * C / lambda3);
	lambda2 = 0.5 * (A - lambda3 + Dlambda21);
	Dlambda32 = lambda3 - lambda2;
	Dlambda31 = Dlambda32 + Dlambda21;

	// ----------------------- //
	// Use Rosetta for Veisq's //
	// ----------------------- //
	// denominators	  
	PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21);
	Xp3 = PiDlambdaInv * Dlambda21;
	Xp2 = -PiDlambdaInv * Dlambda31;

	// numerators
	Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
	Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

	Smm = A - Dmsq21 * Um2sq - Dmsq31 * Um3sq;
	Tmm = Tmm * (1 - Um3sq - Um2sq) + Amatter * (See + Smm - A);

	Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
	Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

	// ------------- //
	// Use NHS for J //
	// ------------- //
	Jmatter = Jmatter * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

	// ----------------------- //
	// Get all elements of Usq //
	// ----------------------- //
	Ue1sq = 1 - Ue3sq - Ue2sq;
	Um1sq = 1 - Um3sq - Um2sq;

	Ut3sq = 1 - Um3sq - Ue3sq;
	Ut2sq = 1 - Um2sq - Ue2sq;
	Ut1sq = 1 - Um1sq - Ue1sq;

	// ----------------------- //
	// Get the kinematic terms //
	// ----------------------- //
	Lover4E = eVsqkm_to_GeV_over4 * L / E;

	D21 = Dlambda21 * Lover4E;
	D32 = Dlambda32 * Lover4E;
	  
	sinD21 = sin(D21);
	sinD31 = sin(D32 + D21);
	sinD32 = sin(D32);

	triple_sin = sinD21 * sinD31 * sinD32;

	sinsqD21_2 = 2 * sinD21 * sinD21;
	sinsqD31_2 = 2 * sinD31 * sinD31;
	sinsqD32_2 = 2 * sinD32 * sinD32;

	// ------------------------------------------------------------------- //
	// Calculate the three necessary probabilities, separating CPC and CPV //
	// ------------------------------------------------------------------- //
	Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
			+ (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
			+ (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
	Pme_CPV = -Jmatter * triple_sin;

	Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
				 + Um3sq * Um1sq * sinsqD31_2
				 + Um3sq * Um2sq * sinsqD32_2);

	Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
				 + Ue3sq * Ue1sq * sinsqD31_2
				 + Ue3sq * Ue2sq * sinsqD32_2);

	// ---------------------------- //
	// Assign all the probabilities //
	// ---------------------------- //
	(*probs_returned)[0][0] = Pee;														// Pee
	(*probs_returned)[0][1] = Pme_CPC - Pme_CPV;										// Pem
	(*probs_returned)[0][2] = 1 - Pee - (*probs_returned)[0][1];  						// Pet

	(*probs_returned)[1][0] = Pme_CPC + Pme_CPV;										// Pme
	(*probs_returned)[1][1] = Pmm;														// Pmm
	(*probs_returned)[1][2] = 1 - (*probs_returned)[1][0] - Pmm;						// Pmt

	(*probs_returned)[2][0] = 1 - Pee - (*probs_returned)[1][0];						// Pte
	(*probs_returned)[2][1] = 1 - (*probs_returned)[0][1] - Pmm;						// Ptm
	(*probs_returned)[2][2] = 1 - (*probs_returned)[0][2] - (*probs_returned)[1][2];	// Ptt
}

void Probability_Vacuum_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, double E, double (*probs_returned)[3][3])
{
	double c13sq, sind, cosd, Jrr, Jvac;
	double Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
	double Lover4E, D21, D31;
	double sinD21, sinD31, sinD32;
	double sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
	double Pme_CPC, Pme_CPV, Pmm, Pee;

	// --------------------------------------------------------------------- //
	// First calculate useful simple functions of the oscillation parameters //
	// --------------------------------------------------------------------- //
	c13sq = 1 - s13sq;

	// Ueisq's
	Ue3sq = s13sq;
	Ue2sq = c13sq * s12sq;

	// Umisq's, Utisq's and Jvac	 
	Um3sq = c13sq * s23sq;
	// Um2sq and Ut2sq are used here as temporary variables, will be properly defined later	 
	Ut2sq = s13sq * s12sq * s23sq;
	Um2sq = (1 - s12sq) * (1 - s23sq);
	  
	Jrr = sqrt(Um2sq * Ut2sq);
	sind = sin(delta);
	cosd = cos(delta);
	Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd;
	Jvac = 8 * Jrr * c13sq * sind;
	
	// ----------------------- //
	// Get all elements of Usq //
	// ----------------------- //
	Ue1sq = 1 - Ue3sq - Ue2sq;
	Um1sq = 1 - Um3sq - Um2sq;

	Ut3sq = 1 - Um3sq - Ue3sq;
	Ut2sq = 1 - Um2sq - Ue2sq;
	Ut1sq = 1 - Um1sq - Ue1sq;

	// ----------------------- //
	// Get the kinematic terms //
	// ----------------------- //
	Lover4E = eVsqkm_to_GeV_over4 * L / E;

	D21 = Dmsq21 * Lover4E;
	D31 = Dmsq31 * Lover4E;
	  
	sinD21 = sin(D21);
	sinD31 = sin(D31);
	sinD32 = sin(D31-D21);

	triple_sin = sinD21 * sinD31 * sinD32;

	sinsqD21_2 = 2 * sinD21 * sinD21;
	sinsqD31_2 = 2 * sinD31 * sinD31;
	sinsqD32_2 = 2 * sinD32 * sinD32;

	// ------------------------------------------------------------------- //
	// Calculate the three necessary probabilities, separating CPC and CPV //
	// ------------------------------------------------------------------- //
	Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
			+ (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
			+ (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
	
	Pme_CPV = -Jvac * triple_sin;

	Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
				 + Um3sq * Um1sq * sinsqD31_2
				 + Um3sq * Um2sq * sinsqD32_2);

	Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
				 + Ue3sq * Ue1sq * sinsqD31_2
				 + Ue3sq * Ue2sq * sinsqD32_2);

	// ---------------------------- //
	// Assign all the probabilities //
	// ---------------------------- //
	(*probs_returned)[0][0] = Pee;														// Pee
	(*probs_returned)[0][1] = Pme_CPC - Pme_CPV;										// Pem
	(*probs_returned)[0][2] = 1 - Pee - (*probs_returned)[0][1];  						// Pet

	(*probs_returned)[1][0] = Pme_CPC + Pme_CPV;										// Pme
	(*probs_returned)[1][1] = Pmm;														// Pmm
	(*probs_returned)[1][2] = 1 - (*probs_returned)[1][0] - Pmm;						// Pmt

	(*probs_returned)[2][0] = 1 - Pee - (*probs_returned)[1][0];						// Pte
	(*probs_returned)[2][1] = 1 - (*probs_returned)[0][1] - Pmm;						// Ptm
	(*probs_returned)[2][2] = 1 - (*probs_returned)[0][2] - (*probs_returned)[1][2];	// Ptt
}

int main()
{
	double L, E, rho, Ye, probs_returned[3][3];
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
	int N_Newton;

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	E = 2.5; // GeV
	rho = 3; // g/cc
	Ye = 0.5;

	// --------------------------------------------------------------------- //
	// Set the number of Newton-Raphson iterations which sets the precision. //
	// 0 is close to the single precision limit and is better than DUNE/HK   //
	// in the high statistics regime. Increasig N_Newton to 1,2,... rapidly  //
	// improves the precision at a modest computational cost                 //
	// --------------------------------------------------------------------- //
	N_Newton = 0;

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = 0.7 * M_PI;
	Dmsq21 = 7.5e-5; // eV^2
	Dmsq31 = 2.5e-3; // eV^2

	// ------------------------------------------ //
	// Calculate all 9 oscillations probabilities //
	// ------------------------------------------ //
	Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);

	// --------------------------- //
	// Print out the probabilities //
	// --------------------------- //
	printf("L = %g E = %g rho = %g\n", L, E, rho);
	printf("Probabilities:\n");
	printf("alpha beta P(nu_alpha -> nu_beta)\n");
	for (int alpha = 0; alpha < 3; alpha++)
	{
		for (int beta = 0; beta < 3; beta++)
		{
			printf("%d %d %g\n", alpha, beta, probs_returned[alpha][beta]);
		} // beta, 3
	} // alpha, 3

}
