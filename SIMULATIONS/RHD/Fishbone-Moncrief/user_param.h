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
#define GRID               UNIFORM
#define PRINT_EVOLV        TRUE

#define Black_Hole_Mass    1.0
#define Black_Hole_Spin    0.0
#define l_0                7.0
#define r_in               30.0

#define MDOT               TRUE
#define MDOT_END           TRUE  
#define MDOT_DATA          TRUE

double rho(double x);
double pre(double x);
double vx1(double x);
double ftheta(double th);
double gtheta(double th);
void Mass_Accretion_Rate(double *B);

double density_0, pressure_0, velocity_0;
double delta;
double cs, Temp;
double rB, tB;
double K_pol;
double Mdot_Mean;

char last[50];
int count, plus, minus, restart_file;
int Mdot_end;
double Mdot_0, Mdot_tprint, Mdot_Max, Mdot_Min, MDOT_ERR, MDOT_TIME;
char eqstate[50];
