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
   
   double r;
   
   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;
   
   polyK = (K-1.)/K*(1.- 1./Rad) ;
      
#if DIM == 1 

   ///////////////////////////
   //-------Bondi-1D--------//
   ///////////////////////////
   
   
   for(i = 0; i <= Nx1; i++)
   {
   	  r = grid.X1[i];
      U(RHO,i) = pow(1. - Rad/r*((r-1.)/(Rad-1.)),1./(K-1.));
      U(PRE,i) = polyK*pow(U(RHO,i),K);
      U(VX1,i) = 0.0;
      
      printf("%e %e \n",r,U(RHO,i));
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
         U(VX1,i,j) = velocity_0;
         U(VX2,i,j) = 0.0;
      }
   }

#elif DIM == 4

   ///////////////////////////
   //------Bondi-2.5D-------//
   ///////////////////////////

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
