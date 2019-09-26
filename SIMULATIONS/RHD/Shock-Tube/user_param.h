#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

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
