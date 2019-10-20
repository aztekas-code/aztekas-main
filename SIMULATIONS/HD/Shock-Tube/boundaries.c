/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 18-10-2019 18:21:45
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
//   for(n = 0; n < eq; n++)
//   {
//      for(i = gc; i <= Nx1-gc; i++)
//      {
//         for(j = gc; j <= Nx2-gc; j++)
//         {
//            if(i > Nx1-gc)
//            {
//               B(n,i,j) = B(n,i-2,j) + ((grid.X1[i] - grid.X1[i-2])/(grid.X1[i-1] - grid.X1[i-2]))*(B(n,i-1,j) - B(n,i-2,j)); 
//            }
//            if(j > Nx2-gc)
//            {
//               B(n,i,j) = B(n,i,j-2) + ((grid.X2[j] - grid.X2[j-2])/(grid.X2[j-1] - grid.X2[j-2]))*(B(n,i,j-1) - B(n,i,j-2)); 
//            }
//         }
//      }
//
//      for(i = Nx1-gc; i >= gc; i--)
//      {
//         for(j = Nx2-gc; j >= gc; j--)
//         {
//            if(i < gc)
//            {
//               B(n,i,j) = B(n,i+2,j) + ((grid.X1[i] - grid.X1[i+2])/(grid.X1[i+1] - grid.X1[i+2]))*(B(n,i+1,j) - B(n,i+2,j)); 
//            }
//            if(j < gc)
//            {
//               B(n,i,j) = B(n,i,j+2) + ((grid.X2[j] - grid.X2[j+2])/(grid.X2[j+1] - grid.X2[j+2]))*(B(n,i,j+1) - B(n,i,j+2)); 
//            }
//         }
//      }
//   }
#endif

   return 0;
}
