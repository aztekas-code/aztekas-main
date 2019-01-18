/* 
 *  aztekas boundaries module
 *  Date of creation: 03-01-2019 02:13:54
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
