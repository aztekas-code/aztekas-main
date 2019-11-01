/*
 * aztekas initial module
 * Date of creation/modification: 27-09-19 11:26:00
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   ///////////////////////////
   //---------Wind----------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) =  density_0;
         U(PRE,i,j) =  pow(U(RHO,i,j),K)/K;
         U(VX1,i,j) =  velocity_0*cos(grid.X2[j]);
         U(VX2,i,j) = -velocity_0*sin(grid.X2[j]);
      }
   }
}
