/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Dust Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:26:36
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
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
         B(VX1,i) = r_dot_0;
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
            B(VX1,i,j) = r_dot_0;
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
            B(VX1,i,j) = r_dot_0;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = phi_dot_0*grid.X1[i]*sin(grid.X2[j]);
         }
      }
   }

#endif

   return 0;
}
