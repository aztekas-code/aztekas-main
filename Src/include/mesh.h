/**
 * @file mesh.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Definitions of the vector, structures and functions needed for the 
 * mesh grid at the center and at the interfaces
 */

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
}grid_;

grid_ grid;

typedef struct
{
   double x[4];
   double lapse;
   double beta_con[3];
   double gamma_con[3][3];
   double dety;
}gauge_;

typedef struct
{
   double dlapse[3];
   double dbeta[3][3];
   double dgam[3][3][3];
}der_gauge_;

double dx1, dx2, dx3;
double dt;
double tmax, cou;

int Mesh(); 
void Surface_Volume();

double TimeStep();
