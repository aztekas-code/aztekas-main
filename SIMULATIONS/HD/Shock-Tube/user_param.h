/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Shock-Tube 
 * Creation Date : 26-09-2019
 * Last Modified : 18-10-2019 17:06:37
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define HORIZONTAL         0
#define VERTICAL           1
#define DIAGONAL           2

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define outflow_x2max      TRUE
#define outflow_x2min      TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM

#define INTERFACE          VERTICAL

double rhol, pl, vx1l, vx2l, vx3l;
double rhor, pr, vx1r, vx2r, vx3r;
double x_0;
