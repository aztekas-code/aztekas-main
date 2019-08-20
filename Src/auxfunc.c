/**
 * @file auxfunc.c
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Helpful functions for \a aztekas.
 */

//Do not erase any of these libraries//
#include"main.h"

int MxV(double *M, double *V, double *L)
{
   int n, m;
   double res=0.0;

   for(m = 0; m < eq; m++)
   {
      for(n = 0; n < eq; n++)
      {
         res += M[m*(eq) + n]*V[n];
      }

      L[m] = res;
      res = 0.0;
   }

   return 0;
}

void RoundGen(double *num)
{
   double r;
   double bla;
   double decnum;
   int expnum;
   
   if(*num != 0.0e+00)
   {
      expnum = (int)floor(log(fabs(*num))/log(10.0));
      decnum = *num/pow(10,(double)expnum);
      decnum = roundf(decnum*1.0e+15)/1.0e+15;
      *num = decnum*pow(10,(double)expnum);
   }
}

void Scalar_Contraction_Range1(double *scalar, double *cov, double *con)
{
   *scalar = cov[0]*con[0] + cov[1]*con[1] + cov[2]*con[2];
}

void Raise_Index_Range1(double *con, double *cov, grid_ *local_grid)
{
   int i, j;

   con[0] = local_grid->gamma_con[0][0]*cov[0] + local_grid->gamma_con[0][1]*cov[1] + local_grid->gamma_con[0][2]*cov[2];
   con[1] = local_grid->gamma_con[1][0]*cov[0] + local_grid->gamma_con[1][1]*cov[1] + local_grid->gamma_con[1][2]*cov[2];
   con[2] = local_grid->gamma_con[2][0]*cov[0] + local_grid->gamma_con[2][1]*cov[1] + local_grid->gamma_con[2][2]*cov[2];


}

void Low_Index_Range1(double *cov, double *con, grid_ *local_grid)
{
   int i, j;

   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      {
         cov[i] += local_grid->gamma_cov[i][j]*con[j];
      }
   }
}

void Low_Index_Range2(double **diag, double **con, grid_ *local_grid)
{
   int i, j, k;

   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      {
         for(k = 0; k < 3; k++)
         {
            diag[i][k] += local_grid->gamma_cov[i][j]*con[j][k];
         }
      }
   }
}

void CheckSimParameters()
{
   printf("\n");
   printf("aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n");
   printf("    a     zz    t    e   e  k  k       a  ss   \n");
   printf("aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n");
   printf("a   a  zz       t    e      k  k   a   a     ss\n");
   printf("aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n");
   printf("\n");
   printf("Running aztekas simulation...\n");
   printf("\n");

   // Print physics used
   if(PHYSICS == HD)  printf("Performing a HD simulation.\n");
   if(PHYSICS == RHD) printf("Performing a RHD simulation.\n");

   // Equation of state
   if(EOS == IDEAL) printf("Ideal equation of state.\n");
   if(EOS == DUST)  printf("Dust.\n");

   // Print MoL-RK order
   printf("Time integration using a second order MoL-Runge Kutta.\n");

   // Print spatial numerical methods, algorithms and parameters.
   // Primitive variable reconstruction.
   if(RECONST == GODUNOV)  printf("Zero-order piecewise reconstruction.\n");
   if(RECONST == MINMOD)   printf("MINMOD reconstruction.\n");
   if(RECONST == MC)       printf("MC reconstruction.\n");
   if(RECONST == SUPERBEE) printf("SUPERBEE reconstruction.\n");
   if(RECONST == WENO5)    printf("WENO5 reconstruction.\n");

   // Flux solver
   if(FLUX == HLL)  printf("HLL Riemann solver.\n");
   if(FLUX == HLLC) printf("HLLC Riemann solver.\n");
}

