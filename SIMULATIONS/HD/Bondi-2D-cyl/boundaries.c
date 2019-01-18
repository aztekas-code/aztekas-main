/* 
 *  aztekas boundaries module
 *  Date of creation: 17-01-2019 20:07:41
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
   REFLECTIVE(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         r = sqrt(pow(X1[i],2.0) + pow(X2[j],2.0));
         if(r <= r_in)
         {
            B[c2(0,i,j)] = 1.0;
            B[c2(1,i,j)] = 1.0;
            B[c2(2,i,j)] = 0.0;
            B[c2(3,i,j)] = 0.0;
            B[c2(4,i,j)] = 0.5;
         }
      }
   }
   return 0;
}
