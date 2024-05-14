#ifndef DMP_H
#define DMP_H

// From http://arxiv.org/abs/1604.08167
// Zeorth order only

void Probability_Matter_LBL_DMP(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, double E, double rho, double Ye, int empty, double (*probs_returned)[3][3]);

#endif
