/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:50:45
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

   PERIODIC(B);
   REFLECTIVE(B);
   for(n = 0; n < eq; n++)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            roundgen(&B[c2(n,i,j)]);
         }
      }
   }

   return 0;
}
