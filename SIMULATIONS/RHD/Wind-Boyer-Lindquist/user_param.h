#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define POLAR              TRUE

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#if POLAR == FALSE
   #define reflective_x2max   TRUE
   #define reflective_x2min   TRUE
#elif POLAR == TRUE
   #define periodic_x2        TRUE
#endif

#define x1min_exc          TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0
#define PRINT_EVOLV        TRUE

#define Black_Hole_Mass    1.0
#define Black_Hole_Spin    0.99
#define Mach               5.0
#define vinf               0.5

double density_0, pressure_0, velocity_0;
