/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Polar Disk
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 17:53:25
 * Created By : Alejandro Aguayo-Ortiz
 */


#include"main.h"

void User_Boundaries(double *B)
{
   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = r_dot_0;
            B(VX2,i,j) = phi_dot_0/grid.X1[i];
         }
      }
   }
}
