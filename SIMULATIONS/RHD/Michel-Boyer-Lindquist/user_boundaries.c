/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Michel Ideal Gas problem
 * Creation Date : 05-12-2019
 * Last Modified : 17-12-2019 13:14:40
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
   int i, j, k, n, cell;
   double r, theta;
   double Delta, Sigma, rho2;
   double M, a;
   double rplus, rminus;

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
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = velocity_0;
            B(3,i,j) = 0.0;
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
}
