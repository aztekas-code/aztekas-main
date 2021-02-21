#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0
#define PRINT_EVOLV        TRUE

#define Black_Hole_Mass    1.0
/**
 * Model 1:    l_0 = 6.657,  c_infty = 0.07071
 * Model 1bis: l_0 = 10.204,  c_infty = 0.07071
 * Model 2:    l_0 = 4.5714, c_infty = 0.02236
 * Model 3:    l_0 = 10.204, c_infty = 0.02236
 * Model 4:    l_0 = 14.286, c_infty = 0.01257
 */
#define l_0                6.657

double rho(double x);
double pre(double x);
double vx1(double x);
double ftheta(double th);
double gtheta(double th);

double density_0, pressure_0, velocity_0;
double delta;
