#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define outflow_x1max 1
#define outflow_x1min 1
#define outflow_x2max 1
#define outflow_x2min 1

#define limiter 'C'
#define riemann 2

#define interface 2
#define logmesh   0

double density_0, pressure_0, velocity_0;
double nl, pl, vx1l, vx2l, vx3l;
double nr, pr, vx1r, vx2r, vx3r;
double x_0;
