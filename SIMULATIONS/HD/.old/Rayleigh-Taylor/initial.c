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

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(X2[j] <= 0.0)
         {
            U[c2(0,i,j)] = 1.0;
            U[c2(1,i,j)] = 2.5 - U[c2(0,i,j)]*X2[j];
            U[c2(2,i,j)] = 0.0;
            U[c2(3,i,j)] = -0.1*(1 + cos(2*M_PI*X1[i]/0.5))*(1 + cos(2*M_PI*X2[j]/1.5))/4;
         }
         if(X2[j] > 0.0)
         {
            U[c2(0,i,j)] = 2.0;
            U[c2(1,i,j)] = 2.5 - U[c2(0,i,j)]*X2[j];
            U[c2(2,i,j)] = 0.0;
            U[c2(3,i,j)] = -0.1*(1 + cos(2*M_PI*X1[i]/0.5))*(1 + cos(2*M_PI*X2[j]/1.5))/4;
         }
      }
   }
}
