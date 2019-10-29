/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 18:20:53
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
#if DIM == 1

   for(int i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
      }
   }

#elif DIM == 2

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0*gtheta(grid.X2[j]);
            B(PRE,i,j) = pow(B(RHO,i,j),K)/K;
         }
      }
   }

#elif DIM == 4

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0*gtheta(grid.X2[j]);
            B(PRE,i,j) = pow(B(RHO,i,j),K)/K;
         }
      }
   }

#endif
}

double gtheta(double th)
{
   double dum = (1.0 - rho_atm)*pow(sin(th),2.0) + rho_atm;
   
   return dum;
}
