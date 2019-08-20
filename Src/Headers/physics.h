#include<mesh.h>

typedef struct
{
	double up[eq+1];
	double um[eq+1];
	double qp[eq+1];
	double qm[eq+1];
	double fp[eq+1];
	double fm[eq+1];
	double lp;
	double lm;
}flx_;

typedef struct
{
	double A[(eq+1)*(eq+1)];
	double S[eq+1];
	double Fp[eq+1];
	double Fm[eq+1];
	double Gp[eq+1];
	double Gm[eq+1];
	double Hp[eq+1];
	double Hm[eq+1];
}vec_;

typedef struct
{
   double e;
   double cs;
   double h;
}eos_;

int funct_A(double *a, double *uu);
int Cons2Prim(double *q, double *u);

void Prim2Cons_All(double *u, double *q);

void Prim2Cons(double *q, double *u, gauge_ local_grid);
void Prim2FluxF(double *f, double *v, double *u, gauge_ local_grid);
void Prim2FluxG(double *f, double *v, double *u, gauge_ local_grid);
void Prim2FluxH(double *f, double *v, double *u, gauge_ local_grid);

void Sources(double *u, vec_ *v, int *I);
void Source_Terms(double *s, double *u, gauge_ local_grid);
void User_Source_Terms(double *s, double *u, gauge_ local_grid);

void EoS(eos_ *e, double *u, gauge_ local_grid);

void Get_Metric_Components(gauge_ *local_grid);
void Scalar_Contraction_Range1(double *scalar, double *cov, double *con);
void Raise_Index_Range1(double *con, double *cov, gauge_ *local_grid);
void Low_Index_Range1(double *cov, double *con, gauge_ *local_grid);
void Low_Index_Range2(double **diag, double **con, gauge_ *local_grid);
