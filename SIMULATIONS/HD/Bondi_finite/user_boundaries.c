/*
 * aztekas boundaries module
 * Date of creation/modification: 26-09-19 11:33:21
 * author: Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
#if DIM == 1

   double r;

   polyK = (K-1.)/K*(1.- 1./Rad) ;
   
   pressure_0 = polyK*pow(density_0,K);

   for(int i = 0; i <= Nx1; i++)
   {

      if(i <= gc)
      {

		  r = grid.X1[i];
	      B(RHO,i) = pow(1. - Rad/r*((r-1.)/(Rad-1.)),1./(K-1.));
    	  B(PRE,i) = polyK*pow(B(RHO,i),K);
    	  B(VX1,i) = 0.0;
	  }

      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
         B(VX1,i) = 0.0;
      }
   }

#elif DIM == 2

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0q;
            B(VX2,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = 0.0;
         }
      }
   }

#endif
}
