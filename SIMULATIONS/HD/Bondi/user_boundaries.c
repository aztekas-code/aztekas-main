/*
 * aztekas boundaries module
 * Date of creation/modification: 26-09-19 11:33:21
 * author: Alejandro Aguayo-Ortiz
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
         B(PRE,i) = 10*pressure_0;
         B(VX1,i) = velocity_0;
      }
   }

#elif DIM == 2

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
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

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
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
}
