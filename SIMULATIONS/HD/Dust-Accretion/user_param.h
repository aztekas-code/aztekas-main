/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Dust-Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:29:13
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

double density_0, pressure_0, r_dot_0, phi_dot_0;
