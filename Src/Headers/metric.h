typedef struct
{
   double gd[4][4];
   double gu[4][4];

   double yd[3][3];
   double yu[3][3];
}metric_;

typedef struct
{
   double lapse;
   double dety;
   double detg;
   double beta[3];
}gauge_;

typedef struct
{
   double dety;
   double gamma[3][3];
}local_metric_;

typedef struct
{
   double *g[DIM+1][DIM+1];
   double *lapse;
   double *beta_u[DIM];
   double *dety;
   double *gamma[DIM][DIM];
   double *beta_d[DIM];
   double *dgamma[DIM][DIM][DIM];
   double *dbeta_u[DIM][DIM];
   double *dlapse[DIM];
}global_metric_;

global_metric_ gmetric;

double dety(double x1, double x2, double x3);
void Metric_Components(local_metric_ *g, double *x);
void Gauge_Components(gauge_ *g, double *x);
