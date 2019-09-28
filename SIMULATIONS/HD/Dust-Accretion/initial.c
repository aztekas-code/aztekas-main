/*
 * File Name : initial.c
 * Description : aztekas initial module for Dust Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:25:19
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   ///////////////////////////
   //-------Bondi-1D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      U(RHO,i) = density_0;
      U(PRE,i) = pressure_0;
      U(VX1,i) = r_dot_0;
   }

#elif DIM == 2

   ///////////////////////////
   //-------Bondi-2D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = r_dot_0;
         U(VX2,i,j) = 0.0;
      }
   }

#elif DIM == 4

   ///////////////////////////
   //-------Bondi-2D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = r_dot_0;
         U(VX2,i,j) = 0.0;
         U(VX3,i,j) = phi_dot_0*grid.X1[i]*sin(grid.X2[j]);
      }
   }

#endif
}
