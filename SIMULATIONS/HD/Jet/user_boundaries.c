/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Jet
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 17:42:46
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
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
}
