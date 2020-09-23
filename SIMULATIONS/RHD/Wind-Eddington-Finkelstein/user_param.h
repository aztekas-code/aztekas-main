#include"macros.h"

// Define boundaries TRUE = 1 , FALSE = 0
// Set as cond_X whereas
// cond : {outflow,reflective,periodic,inflow}
// X : {x1max,x1min,x2max,x2min,x3max,x3min}

#define POLAR              FALSE

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
#define GRID               UNIFORM

#define Black_Hole_Mass    1.0
#define Black_Hole_Spin    0.0

#define MDOT               TRUE
#define MDOT_END           FALSE
#define MDOT_ERR           1.0e-4 // only if MDOT_END TRUE
#define MDOT_TIME          1.0
int Mdot_end;
double Mdot_0;
double Mdot_tprint;
double Mdot_Mean;
double Mdot_Max;
double Mdot_Min;
char last[50];
int count, plus, minus;
int restart_file;

void Mass_Accretion_Rate(double *B);
double density_inf, pressure_inf, velocity_inf, Mach;
