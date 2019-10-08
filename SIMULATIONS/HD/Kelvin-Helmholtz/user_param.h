/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:08:27
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define periodic_x1        TRUE
#define periodic_x2        TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM

double rhod, pd, vx1d, vx2d, vx3d;
double rhou, pu, vx1u, vx2u, vx3u;
double x_0;
double eta;
