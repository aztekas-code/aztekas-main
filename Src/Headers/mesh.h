typedef struct
{
   double time;
   double *X1;
   double *X1p, *X1m;
   double *X2;
   double *X2p, *X2m;
   double *X3;
   double *X3p, *X3m;

   double *S1p, *S1m;
   double *S2p, *S2m;
   double *S3p, *S3m;

   double x[4];

   double lapse;
   double *beta_cov;
   double *beta_conv;
   double *gamma_cov;
   double *gamma_con;

   double *g_cov;
   double *g_con;
}grid_;

double dx1, dx2, dx3;
double dt;
double tmax, timefile, cou;

grid_ grid;

int Mesh(); 
void Surface_Volume();

double TimeStep();
