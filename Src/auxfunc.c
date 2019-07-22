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

   // Print MoL-RK order
   printf("Time integration using a second order MoL-Runge Kutta.\n");

   // Print spatial numerical methods, algorithms and parameters.
   // Primitive variable reconstruction.
   if(RECONST == GODUNOV)  printf("Using a zero-order piecewise reconstruction.\n");
   if(RECONST == MINMOD)   printf("Using a MINMOD reconstruction.\n");
   if(RECONST == MC)       printf("Using a MC reconstruction.\n");
   if(RECONST == SUPERBEE) printf("Using a SUPERBEE reconstruction.\n");
   if(RECONST == WENO5)    printf("Using a WENO5 reconstruction.\n");
}
