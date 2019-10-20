#include<mesh.h>

typedef struct
{
	double U[eq+1];
	double Q[eq+1];
   double F[eq+1];
   double L[eq+1];
	double A[(eq+1)*(eq+1)];
	double Fp[eq+1];
	double Fm[eq+1];
	double Gp[eq+1];
	double Gm[eq+1];
	double Hp[eq+1];
	double Hm[eq+1];
	double S[eq+1];
}rhs_;

typedef struct
{
   double e;
   double cs;
   double h;
   double dhdrho;
   double dhdp;
}eos_;

int funct_A(double *a, double *uu);
int Cons2Prim(double *q, double *u);

void Prim2Cons_All(double *u, double *q);

void Prim2Cons(double *q, double *u, gauge_ *local_grid);
void Prim2FluxF(double *f, double *v, double *u, gauge_ *local_grid);
void Prim2FluxG(double *f, double *v, double *u, gauge_ *local_grid);
void Prim2FluxH(double *f, double *v, double *u, gauge_ *local_grid);

void Prim2Sources(double *s, int *I);
void Source_Terms(double *s, double *u, gauge_ *local_grid);
void User_Source_Terms(double *s, double *u, gauge_ *local_grid);

void Matrix_A(double *a, double *u, gauge_ local_grid);

void EoS(eos_ *e, double *u, gauge_ *local_grid);

void Get_Metric_Components(gauge_ *local_grid);
void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid);
void Scalar_Contraction_Range1(double *scalar, double *cov, double *con);
void Raise_Index_Range1(double *con, double *cov, gauge_ *local_grid);
void Low_Index_Range1(double *cov, double *con, gauge_ *local_grid);
void Low_Index_Range2(double **diag, double **con, gauge_ *local_grid);
