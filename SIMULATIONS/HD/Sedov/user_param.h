/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Sedov Blast Wave
 * Creation Date : 26-09-2019
 * Last Modified : 15-10-2019 17:22:08
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define GRAVITY            NONE
#define USER_SOURCE_TERMS  FALSE

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min      TRUE

#define RECONST            MC
#define FLUX               HLLC
#define GRID               UNIFORM

double rho_0, p_0, E_0;
double x_0;
