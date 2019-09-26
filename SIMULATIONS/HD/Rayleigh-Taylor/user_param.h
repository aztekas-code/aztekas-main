/*
 * aztekas user parameters header file
 * Date of creation/modification: 25-09-19 23:58:32
 * author: Alejandro Aguayo-Ortiz
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
