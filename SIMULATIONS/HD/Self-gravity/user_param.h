/*
 * aztekas user parameters header file
 * Date of creation/modification: 26-09-19 11:29:53
 * author: Alejandro Aguayo-Ortiz
 */

//#include"macros.h"

#define outflow_x1max      FALSE
#define outflow_x1min      FALSE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MINMOD
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0
#define PRINT_EVOLV        TRUE
#define selfgravity        TRUE


double density_0, pressure_0, velocity_0, Rad, polyK, r_acc;

// mass that has been accreted (r<r_acc)
double mass_acc = 0.0;
// Total gas cloud's mass
double mass_tot;
double *mass;
