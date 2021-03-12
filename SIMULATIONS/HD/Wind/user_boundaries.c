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
            B(RHO,i,j) =  density_0;
            B(PRE,i,j) =  pressure_0;//pow(B(RHO,i,j),K)/K;
            B(VX1,i,j) =  velocity_0*cos(grid.X2[j]);
            B(VX2,i,j) = -velocity_0*sin(grid.X2[j]);
         }

         if(i <= gc)
         {
            B(RHO,i,j) = B(RHO,gc+1,j);
            B(PRE,i,j) = B(PRE,gc+1,j);
            B(VX1,i,j) = B(VX1,gc+1,j);
            B(VX2,i,j) = B(VX2,gc+1,j);
         }
      }
   }
}
