#include"macros.h"

#define POLAR              TRUE

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#if POLAR == FALSE
   #define reflective_x2max   TRUE
   #define reflective_x2min   TRUE
#elif POLAR == TRUE
   #define periodic_x2        TRUE
#endif

#define RECONST            MC
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0

#define Black_Hole_Mass    1.0
#define Black_Hole_Spin    0.995
#define Mach               5.0
#define vinf               0.5

double density_0, pressure_0, velocity_0;
