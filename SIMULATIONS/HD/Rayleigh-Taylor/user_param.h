/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Shock-Tube
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:08:50
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define periodic_x1        TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM

double rhod, pd, vx1d, vx2d, vx3d;
double rhou, pu, vx1u, vx2u, vx3u;
double x_0;
double eta;
