/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 26-09-2019 23:57:13
 * Created By :
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

#if DIM == 1
   Outflow(B);

#elif DIM == 2
   Outflow(B);
   
   // Linear extrapolation. Important for the diagonal shock tube
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(n = 0; n < eq; n++)
         {
            if(i > Nx1-gc)
            {
               B(n,i,j) = B(n,i-2,j) + ((grid.X1[i] - grid.X1[i-2])/(grid.X1[i-1] - grid.X1[i-2]))*(B(n,i-1,j) - B(n,i-2,j)); 
            }
            if(j > Nx2-gc)
            {
               B(n,i,j) = B(n,i,j-2) + ((grid.X2[j] - grid.X2[j-2])/(grid.X2[j-1] - grid.X2[j-2]))*(B(n,i,j-1) - B(n,i,j-2)); 
            }
         }
      }
   }

   for(i = Nx1; i >= 0; i--)
   {
      for(j = Nx2; j <= 0; j--)
      {
         for(n = 0; n < eq; n++)
         {
            if(i > Nx1-gc)
            {
               B(n,i,j) = B(n,i+2,j) + ((grid.X1[i] - grid.X1[i+2])/(grid.X1[i+1] - grid.X1[i+2]))*(B(n,i+1,j) - B(n,i+2,j)); 
            }
            if(j > Nx2-gc)
            {
               B(n,i,j) = B(n,i,j+2) + ((grid.X2[j] - grid.X2[j+2])/(grid.X2[j+1] - grid.X2[j+2]))*(B(n,i,j+1) - B(n,i,j+2)); 
            }
         }
      }
   }
   
#endif

   return 0;
}
