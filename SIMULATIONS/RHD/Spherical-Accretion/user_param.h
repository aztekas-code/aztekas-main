/* 
 *  aztekas header paremeters module
 *  Date of creation: 27-11-2019 10:44:40
 *  author: Alejandro Aguayo-Ortiz 
 */

#define outflow_x1max           TRUE
#define outflow_x1min           TRUE
#define reflective_x2max        TRUE
#define reflective_x2min        TRUE

#define RECONST                 MC
#define FLUX                    HLL
#define GRID                    LOGMESH
#define PRINT_EVOLV             TRUE
#define lfac                    1.0

//#define x1min_exc               TRUE
//#define x1max_exc               TRUE

#define Black_Hole_Mass         1.0

#define MDOT                    TRUE
#define MDOT_END                TRUE
#define MDOT_DATA               1.0e+04
int Mdot_end;
double Mdot_0;
double Mdot_tprint;

void Mass_Accretion_Rate(double *B);
double rs(double K, double c_s);

double Black_Hole_Spin;
double Mdot_Mean;
double cs, Temp;
double rB, tB;
double density_0, velocity_0;
char last[50];

int count;
int plus;
int minus;
int restart_file;
double Mdot_Max;
double Mdot_Min;
double MDOT_ERR;
double MDOT_TIME;
char eqstate[50];
