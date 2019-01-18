#include"const.h"
#include"cond.h"
#include"limiters.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

double *U, *U1, *U2, *U3;
double *Q, *Q1, *Q2, *Q3;
double *X1;
double *X1p, *X1m;
double *X2;
double *X2p, *X2m;
double *X3;
double *X3p, *X3m;

double x1, x2, x3;
double dx1, dx2, dx3;
double dt, time;
double tmax, timefile, cou;

double start, delta;
double K;

//Mesh
int Nx1, Nx2, Nx3;

double x1max, x2max, x3max;
double x1min, x2min, x3min;

//RIEMANN
double nl, pl, vx1l, vx2l, vx3l;
double nr, pr, vx1r, vx2r, vx3r;
double x_0;

//JET
double r_jet;
double z_jet;

double n_jet, p_jet, vx1_jet, vx2_jet, vx3_jet;
double n_atm, p_atm, vx1_atm, vx2_atm, vx3_atm;

//Spherical Accretion
int rho_boundary;
double r_out, r_in;
double theta_0, delta_theta;
double density_0, pressure_0, velocity_0;

//Paramfile
char paramfile_name[50], outputdirectory[50], outputfile[50];
char restartfile[50];
int read_parameters_file(char const *paramfile_name);
int restart_simulation, restart_filecount;

typedef struct
{
	double ux1p[2*eq];
	double ux1m[2*eq];
	double sx1[2*eq];
	double ux2p[2*eq];
	double ux2m[2*eq];
	double sx2[2*eq];
	double ux3p[2*eq];
	double ux3m[2*eq];
	double sx3[2*eq];
	double ux[2*eq];
}lim_;

typedef struct
{
	double up[eq+1];
	double um[eq+1];
	double qp[eq+1];
	double qm[eq+1];
	double fp[eq+1];
	double fm[eq+1];
	double gp[eq+1];
	double gm[eq+1];
	double hp[eq+1];
	double hm[eq+1];
	double lp;
	double lm;
}flx_;

typedef struct
{
	double A[(eq+1)*(eq+1)];
   double Q[eq+1];
   double Q1[eq+1];
   double Q2[eq+1];
	double S[eq+1];
   double F[eq+1];
   double G[eq+1];
   double H[eq+1];
	double Fp[eq+1];
	double Fm[eq+1];
	double Gp[eq+1];
	double Gm[eq+1];
	double Hp[eq+1];
	double Hm[eq+1];
}vec_;

double LIMITER(double A, double B, char r);

double GODUNOV(double A, double B);

double MAXMOD(double A, double B);

double MINMOD(double A, double B);

double MC(double A, double B);

double SUPERBEE(double A, double B);

double WENO5(double v1, double v2, double v3, double v4, double v5);

void allocateArray();

void new_SIZE();

int MESH(); 

void INITIAL();

double TIMESTEP();

int BOUNDARIES(double *B);

int PrintValues(double *tprint, double *dtprint, int *itprint);

int Output1(int *itprint);

int Output2(int *itprint);

int Output3(int *itprint);

int INTEGRATION();

int RK1D(double *u, double *q, double *q1, double *q2, int order);

int RK2D(double *u, double *q, double *q1, double *q2, int order);

int RK3D(double *u, double *q, double *q1, double *q2, int order);

int FLUX1D(vec_ *v, lim_ *l, int *I);
                                   
int FLUX2D(vec_ *v, lim_ *l, int *I);
                                   
int FLUX3D(vec_ *v, lim_ *l, int *I);

int HLL(double *F, flx_ *f, int x);

int HLLC(double *F, flx_ *f, int x);

int AMATRIX1D(double *u, vec_ *v, int *I);
                                                            
int AMATRIX2D(double *u, vec_ *v, int *I);
                                              
int AMATRIX3D(double *u, vec_ *v, int *I);

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I);

int c1(int n, int i);

int c2(int n, int i, int j);

int c3(int n, int i, int j, int k);

int RECONST1D(double *u, char r, lim_ *l, int *I);

int RECONST2D(double *u, char r, lim_ *l, int *I);

int RECONST3D(double *u, char r, lim_ *l, int *I);

int MxV(double *M, double *V, double *L);

void roundgen(double *num);
