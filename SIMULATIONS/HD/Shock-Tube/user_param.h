/*
 * aztekas user parameters header file
 * Date of creation/modification: 25-09-19 22:46:36
 * author: Alejandro Aguayo-Ortiz
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

#define INTERFACE          HORIZONTAL

double rhol, pl, vx1l, vx2l, vx3l;
double rhor, pr, vx1r, vx2r, vx3r;
double x_0;
