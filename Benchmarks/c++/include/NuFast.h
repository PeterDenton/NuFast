#ifndef NuFast_H
#define NuFast_H

void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, double E, double rho, double Ye, int N_Newton, double (*probs_returned)[3][3]);
void Probability_Vacuum_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, double E, double empty1, double empty2, int empty3, double (*probs_returned)[3][3]);

#endif
