/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 15-10-2019 20:21:33
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

   Outflow(B);
   Reflection(B);

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(0,i) = density_0;
         B(1,i) = pressure_0;
         B(2,i) = velocity_0;
      }
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B(0,i,j) = density_0*gtheta(grid.X2[j]);
            B(1,i,j) = pow(B(0,i,j),K)/K;
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
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = velocity_0;
            B(3,i,j) = 0.0;
            B(4,i,j) = 0.0;
         }
      }
   }

#endif

   return 0;
}

double gtheta(double th)
{
   double dum = (1.0 - rho_atm)*pow(sin(th),2.0) + rho_atm;
   
   return dum;
}
