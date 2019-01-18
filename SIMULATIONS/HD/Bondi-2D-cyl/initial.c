/* 
 *  aztekas initial module
 *  Date of creation: 17-01-2019 20:07:41
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
   double r;

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         r = sqrt(pow(X1[i],2.0) + pow(X2[j],2.0));
         U[c2(0,i,j)] = 1.0;
         U[c2(1,i,j)] = 1.0;
         U[c2(2,i,j)] = 0.0;
         U[c2(3,i,j)] = 0.0;
         U[c2(4,i,j)] = 1.0/r;
      }
   }
}
