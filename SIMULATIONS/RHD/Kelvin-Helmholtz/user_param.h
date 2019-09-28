/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Relativistic Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:47:25
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define periodic_x1        TRUE
#define periodic_x2        TRUE

#define RECONST            WENO5
#define FLUX               HLL
#define GRID               UNIFORM

double rhod, pd, vx1d, vx2d, vx3d;
double rhou, pu, vx1u, vx2u, vx3u;
double x_0;
