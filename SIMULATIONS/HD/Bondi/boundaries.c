/* 
 *  aztekas boundaries module
 *  Date of creation: 17-01-2019 16:10:43
 *  author: Alejandro Aguayo Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

int BOUNDARIES(double *B)
{
   int n, i, j, k, cell;
   double r;

   OUTFLOW(B);

#if dim == 1

   for(i = 0; i <= Nx1; i++)
   {
      r = X1[i];

      if(i == Nx1-gc)
      {
         B[c1(0,i)] = density_0;
         B[c1(1,i)] = pressure_0;
         B[c1(2,i)] = velocity_0;
      }
      if(i > Nx1-gc)
      {
         for(n = 0; n < eq; n++)
         {
            B[c1(n,i)] = B[c1(n,i-2)] + (B[c1(n,i-1)] - B[c1(n,i-2)])*(X1[i] - X1[i-2])/(X1[i-1] - X1[i-2]);
         }
      }
      if(i < gc)
      {
         for(n = 0; n < eq; n++)
         {
            B[c1(n,i)] = B[c1(n,i+2)] + (B[c1(n,i+1)] - B[c1(n,i+2)])*(X1[i] - X1[i+2])/(X1[i+1] - X1[i+2]);
         }
      }
   }

   for(i = 0; i <= Nx1; i++)
   {
      for(n = 0; n < eq; n++)
      {
         roundgen(&B[c1(n,i)]);
      }
   }

#endif

   return 0;
}
