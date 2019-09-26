/*
 * aztekas initial module
 * Date of creation: 26-09-2019 00:09:12
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   for(i = 0; i <= Nx1; i++)
   {
      if(fabs(grid.X1[i]) < x_0*dx1)
      {
         U(RHO,i) = rho_0;
         U(PRE,i) = 3*(K-1)*E_0/(4*M_PI*rho_0*pow((x_0)*dx1,3.0));
         U(VX1,i) = 0.0;
      }
      else
      {
         U(RHO,i) = rho_0;
         U(PRE,i) = p_0;
         U(VX1,i) = 0.0;
      }
   }

#elif DIM == 2
   
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(fabs(grid.X1[i]) < x_0*dx1)
         {
            U(RHO,i,j) = rho_0;
            U(PRE,i,j) = 3*(K-1)*E_0/(4*M_PI*rho_0*pow(x_0*dx1,3.0));
            U(VX1,i,j) = 0.0;
            U(VX2,i,j) = 0.0;
         }
         else
         {
            U(RHO,i,j) = rho_0;
            U(PRE,i,j) = p_0;
            U(VX1,i,j) = 0.0;
            U(VX2,i,j) = 0.0;
         }
      }
   }
#endif
}
