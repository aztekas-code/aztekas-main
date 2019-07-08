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

typedef struct
{
   double e;
   double cs;
}eos_;

int funct_A(double *a, double *uu);
int funct_Q2U(double *q, double *u);

int GAUGE(double *a, double g1, double g2, double g3);

void Prim2Cons(double *a, double *uu);
void Prim2Cons_All(double *u, double *q);
void Prim2FluxF(double *f, double *v, double *uu);
void Prim2FluxG(double *f, double *v, double *uu);
void Prim2FluxH(double *f, double *v, double *uu);
int funct_Dm(double *a, double *uu);
int funct_Dn(double *a, double *uu);
int funct_Do(double *a, double *uu);
void Sources(double *u, vec_ *v, int *I);
void Source_Terms(double *s, double *uu, double *x);

void EoS_Ideal(eos_ *e, double *u, double *x);
