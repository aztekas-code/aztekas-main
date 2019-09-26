/*
 * azetkas user parameters header file
 * Date of creation/modification: 26-09-19 00:22:25
 * author: Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define reflective_x1min   TRUE
#define outflow_x2max      TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM

double rho_jet, p_jet, vx1_jet, vx2_jet, vx3_jet;
double rho_atm, p_atm, vx1_atm, vx2_atm, vx3_atm;
double r_jet, z_jet;
