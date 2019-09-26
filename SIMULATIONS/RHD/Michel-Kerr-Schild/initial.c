/*
 * aztekas initial module
 * Date of creation/modification: 26-09-19 11:56:23
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

#if DIM == 1 

   ////////////////////////////
   //-------Michel-1D--------//
   ////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      U(RHO,i) = density_0;
      U(PRE,i) = pressure_0;
      U(VX1,i) = velocity_0;
   }

#elif DIM == 2

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = velocity_0;
         U(VX2,i,j) = 0.0;
      }
   }

#elif DIM == 4

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = velocity_0;
         U(VX2,i,j) = 0.0;
         U(VX3,i,j) = 0.0;
      }
   }

#endif
}
