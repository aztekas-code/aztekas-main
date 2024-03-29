/*
 * aztekas boundaries module
 * Date of creation/modification: 26-09-19 11:33:21
 * author: Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
         {
            B(RHO,i,j) =  density_inf;
            B(PRE,i,j) =  pressure_inf;
            B(VX1,i,j) =  velocity_inf*cos(grid.X2[j]);
            B(VX2,i,j) = -velocity_inf*sin(grid.X2[j]);
         }

         if(i <= gc)
         {
            B(RHO,i,j) = B(RHO,gc+1,j);
            B(PRE,i,j) = B(PRE,gc+1,j);
            B(VX1,i,j) = -fabs(B(VX1,gc+1,j));
            B(VX2,i,j) = B(VX2,gc+1,j);
         }

         //if(B(RHO,i,j) < density_inf)
         //{
         //   B(RHO,i,j) = density_inf;
         //}
      }
   }
}
