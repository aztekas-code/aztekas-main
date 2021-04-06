/*
 * aztekas user parameters header file
 * Date of creation/modification: 26-09-19 11:29:53
 * author: Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MINMOD
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0
#define PRINT_EVOLV        TRUE
#define HELMHOLTZ_COMP     CO1
#define ANALYSIS           TRUE

double density_inf, pressure_inf, velocity_inf, temperature_inf, Mach_inf;
double dens_units, vel_units, temp_units;
