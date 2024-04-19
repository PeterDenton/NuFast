/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#ifndef GLB_OSZPROB_H
#define GLB_OSZPROB_H 1

#include <complex.h>

typedef struct osc_params_angles_type {double theta12, theta13, theta23, delta, dms, dma;} osc_params_angles;

/* User-defined glb_probability_matrix and glb_set/get_oscillation_parameters */
typedef int (*glb_probability_matrix_function)(double P[3][3], int cp_sign, double E,
                  const double *length, const double *density);
typedef int (*glb_set_oscillation_parameters_function)(osc_params_angles p);
typedef int (*glb_get_oscillation_parameters_function)(osc_params_angles p);

#define GLB_V_FACTOR 7.63247e-14   /* Conversion factor for matter potentials */

#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
                                                      /* which vacuum algorithms are used   */
#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */

/* Macros */
#define SQR(x)      ((x)*(x))                        /* x^2   */
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */

/* Functions */
int glb_init_probability_engine();
int glb_free_probability_engine();
int glb_set_oscillation_parameters(osc_params_angles opa);
int glb_get_oscillation_parameters(osc_params_angles opa);
int glb_hamiltonian_cd(double E, double V, int cp_sign);
int glb_S_matrix_cd(double E, double L, double V, int cp_sign);
int glb_probability_matrix(double P[3][3], int cp_sign, double E,
    const double *length, const double *density);

#endif /* GLB_OSZPROB_H */


