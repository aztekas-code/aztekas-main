/*
 * aztekas boundaries module for ideal gas accretion (Michel)
 * Date of creation/modification: 26-09-19 12:26:02
 * author: Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;
   double r, theta;
   double Delta, Sigma, rho2;
   double M, a;
   double rplus, rminus;
   double alpha, betar, grr, grp, gpp;
   double Vr, Vp, vr, vp;

   Outflow(B);
#if POLAR == FALSE
   Reflection(B);
#elif POLAR == TRUE
   Periodic(B);
#endif

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
         B(VX1,i) = velocity_0;
      }
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {

            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = 0.0;
         }
      }
   }

#endif

   return 0;
}
