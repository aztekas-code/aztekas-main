/*
 * aztekas user parameters header file
 * Date of creation/modification: 26-09-19 11:29:53
 * author: Alejandro Aguayo-Ortiz
 */

//#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLLC
#define GRID               UNIFORM

double density_0, pressure_0, velocity_0;
