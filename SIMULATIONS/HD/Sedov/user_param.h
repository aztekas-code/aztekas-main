/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Sedov Blast Wave
 * Creation Date : 26-09-2019
 * Last Modified : 27-09-2019 00:08:46
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min      TRUE

#define RECONST            MC
#define FLUX               HLLC
#define GRID               UNIFORM

double rho_0, p_0, E_0;
double x_0;
