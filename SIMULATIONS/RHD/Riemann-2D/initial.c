/*
 * File Name : initial.c
 * Description : aztekas initial module for Relativistic Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:56:15
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

   ////////////////
   // RIEMANN-2D //
   ////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         // TOP-LEFT
         if(fabs(grid.X2[j]) > x_0 && fabs(grid.X1[i]) <= x_0)
         {
            U(RHO,i,j) = rhotl;
            U(PRE,i,j) = ptl;
            U(VX1,i,j) = vx1tl;
            U(VX2,i,j) = vx2tl;
         }
         // TOP-RIGHT
         else if(fabs(grid.X2[j]) > x_0 && fabs(grid.X1[i]) >  x_0)
         {
            U(RHO,i,j) = rhotr;
            U(PRE,i,j) = ptr;
            U(VX1,i,j) = vx1tr;
            U(VX2,i,j) = vx2tr;
         }
         // BOTTOM-LEFT
         else if(fabs(grid.X2[j]) <=  x_0 && fabs(grid.X1[i]) <= x_0)
         {
            U(RHO,i,j) = rhobl;
            U(PRE,i,j) = pbl;
            U(VX1,i,j) = vx1bl;
            U(VX2,i,j) = vx2bl;
         }
         // BOTTOM-RIGHT
         else if(fabs(grid.X2[j]) <=  x_0 && fabs(grid.X1[i]) > x_0)
         {
            U(RHO,i,j) = rhobr;
            U(PRE,i,j) = pbr;
            U(VX1,i,j) = vx1br;
            U(VX2,i,j) = vx2br;
         }
      }
   }
   //////////////////////////////////
}
