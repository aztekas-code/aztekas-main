/*
 * File Name : initial.c
 * Description : aztekas initial module for Rayleigh-Taylor
 * Creation Date : 27-09-2019
 * Last Modified : 07-10-2019 21:50:18
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   int n, i, j, k, cell;
   double Lx1, Lx2;

   Lx1 = fabs(x1max - x1min);
   Lx2 = fabs(x2max - x2min);

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X2[j] <= 0.0)
         {
            U(RHO,i,j) =  rhod;
            U(PRE,i,j) =  p_0 - U(RHO,i,j)*g*grid.X2[j];
            U(VX1,i,j) =  vx1d;
            U(VX2,i,j) = -eta*(1.0 + cos(2.0*M_PI*grid.X1[i]/Lx1))*(1 + cos(2.0*M_PI*grid.X2[j]/Lx2))/4.0;
         }
         if(grid.X2[j] > 0.0)
         {
            U(RHO,i,j) =  rhou;
            U(PRE,i,j) =  p_0 - U(RHO,i,j)*g*grid.X2[j];
            U(VX1,i,j) =  vx1u;
            U(VX2,i,j) = -eta*(1.0 + cos(2.0*M_PI*grid.X1[i]/Lx1))*(1.0 + cos(2.0*M_PI*grid.X2[j]/Lx2))/4.0;
         }
      }
   }
}
