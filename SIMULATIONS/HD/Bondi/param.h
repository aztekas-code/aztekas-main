#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define outflow_x1max 1
#define outflow_x1min 1
#define reflective_x2max 1
#define reflective_x2min 1

#define limiter 'C'
#define riemann 1

#define logmesh  0

double density_0, pressure_0, velocity_0;
