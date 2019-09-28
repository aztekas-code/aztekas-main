/*
 * File Name : initial.c
 * Description : aztekas initial module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 26-09-2019 23:55:50
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

   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      if(grid.X1[i] < x_0)
      {
         U(RHO,i) = rhol;
         U(PRE,i) = pl;
         U(VX1,i) = vx1l;
      }
      else
      {
         U(RHO,i) = rhor;
         U(PRE,i) = pr;
         U(VX1,i) = vx1r;
      }
   }

   printf("%e\n",K);

#elif DIM == 2

   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////

   #if INTERFACE == HORIZONTAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] < x_0)
         {
            U(RHO,i,j) = rhol;
            U(PRE,i,j) = pl;
            U(VX1,i,j) = vx1l;
            U(VX2,i,j) = vx2l;
         }
         else 
         {
            U(RHO,i,j) = rhor;
            U(PRE,i,j) = pr;
            U(VX1,i,j) = vx1r;
            U(VX2,i,j) = vx2r;
         }
      }
   }
   #elif INTERFACE == VERTICAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X2[j] > x_0)
         {
            U(RHO,i,j) = rhol;
            U(PRE,i,j) = pl;
            U(VX1,i,j) = vx1l;
            U(VX2,i,j) = vx2l;
         }
         else 
         {
            U(RHO,i,j) = rhor;
            U(PRE,i,j) = pr;
            U(VX1,i,j) = vx1r;
            U(VX2,i,j) = vx2r;
         }
      }
   }
   #elif INTERFACE == DIAGONAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] + grid.X2[j] - 1 < 0.0)
         {
            U(RHO,i,j) = rhol;
            U(PRE,i,j) = pl;
            U(VX1,i,j) = vx1l;
            U(VX2,i,j) = vx2l;
         }
         else 
         {
            U(RHO,i,j) = rhor;
            U(PRE,i,j) = pr;
            U(VX1,i,j) = vx1r;
            U(VX2,i,j) = vx2r;
         }
      }
   }
   #endif
#endif
}
