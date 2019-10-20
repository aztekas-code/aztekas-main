/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Jet
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:10:45
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int n, i, j, k, cell;

   Outflow(B);
   Reflection(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] <= r_jet && grid.X2[j] <= z_jet)
         {
            B(RHO,i,j) = rho_jet;
            B(PRE,i,j) = p_jet;
            B(VX1,i,j) = vx1_jet;
            B(VX2,i,j) = vx2_jet;
         }
      }
   }
   

   return 0;
}
