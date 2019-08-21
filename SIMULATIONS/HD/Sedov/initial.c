/*
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
         U(0,i) = n_0;
         U(1,i) = 3*(K-1)*E_0/(4*M_PI*n_0*pow((x_0)*dx1,3.0));
         U(2,i) = 0.0;
      }
      else
      {
         U(0,i) = n_0;
         U(1,i) = p_0;
         U(2,i) = 0.0;
      }
   }

#elif DIM == 2
   
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(fabs(grid.X1[i]) < x_0*dx1)
         {
            U(0,i,j) = n_0;
            U(1,i,j) = 3*(K-1)*E_0/(4*M_PI*n_0*pow(x_0*dx1,3.0));
            U(2,i,j) = 0.0;
            U(3,i,j) = 0.0;
         }
         else
         {
            U(0,i,j) = n_0;
            U(1,i,j) = p_0;
            U(2,i,j) = 0.0;
            U(3,i,j) = 0.0;
         }
      }
   }
#endif
}
