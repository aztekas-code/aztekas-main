/*
 * File Name : initial.c
 * Description : aztekas initial module for Relativistic Jet
 * Creation Date : 27-09-2019
 * Last Modified : 07-10-2019 17:00:55
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] <= r_jet && grid.X2[j] <= z_jet)
         {
            U(RHO,i,j) = rho_jet;
            U(PRE,i,j) = p_jet;
            U(VX1,i,j) = vx1_jet;
            U(VX2,i,j) = vx2_jet;
         }
         else
         {
            U(RHO,i,j) = rho_atm;
            U(PRE,i,j) = p_atm;
            U(VX1,i,j) = vx1_atm;
            U(VX2,i,j) = vx2_atm;
         }
      }
   }
}
