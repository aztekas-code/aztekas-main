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
   double beta_cov[3];
   double beta_con[3];
   double gamma_cov[3][3];
   double gamma_con[3][3];
   double dety;

   double g_cov[4][4];
   double g_con[4][4];
}grid_;

double dx1, dx2, dx3;
double dt;
double tmax, timefile, cou;

grid_ grid;

int Mesh(); 
void Surface_Volume();

double TimeStep();
