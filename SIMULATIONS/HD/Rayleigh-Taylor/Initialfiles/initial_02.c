/* 
 *  aztekas initial module
 *  Date of creation: 02-01-2019 12:50:45
 *  author: Alejandro Aguayo Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include"main.h"
#include"param.h"

void INITIAL(double *dtprint)
{
   int n, i, j, k, cell;
   double max, min, range, div;

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         max = 1.0;
         min = -1.0;
         range = max - min;
         div = RAND_MAX/range;
         if(X2[j] <= 0.0)
         {
            U[c2(0,i,j)] = 1.0;
            U[c2(1,i,j)] = 2.5 - U[c2(0,i,j)]*X2[j];
            U[c2(2,i,j)] = 0.0;
            U[c2(3,i,j)] = 0.01*(min + rand()/div);
         }
         if(X2[j] > 0.0)
         {
            U[c2(0,i,j)] = 2.0;
            U[c2(1,i,j)] = 2.5 - U[c2(0,i,j)]*X2[j];
            U[c2(2,i,j)] = 0.0;
            U[c2(3,i,j)] = 0.01*(min + rand()/div);
         }
      }
   }
}
