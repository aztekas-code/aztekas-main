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
#define ANALYSIS           FALSE

#define MDOT               TRUE
#define MDOT_END           FALSE
#define MDOT_TIME          1.0
#define MDOT_ERR           1.0e-04

double density_inf, pressure_inf, velocity_inf, temperature_inf, Mach_inf;
double dens_units, vel_units, temp_units;

double rho_atm;
char last[50];
int count, plus, minus, restart_file;
int Mdot_end;
double Mdot_Max, Mdot_Min;
double Mdot_0, Mdot_tprint, Mdot_Mean;

void Mass_Accretion_Rate(double *B);
double gtheta(double th);
