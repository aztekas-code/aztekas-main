/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 09:58:40
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               LOGMESH


double density_0, pressure_0, velocity_0;
double rho_atm;

double gtheta(double th);
