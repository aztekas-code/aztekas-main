/* 
 *  aztekas initial module
 *  Date of creation: 03-01-2019 02:13:54
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
         if(X1[i] + X2[j] < x_0)
         {
            U[c2(0,i,j)] = nl;
            U[c2(1,i,j)] = pl;
            U[c2(2,i,j)] = vx1l;
            U[c2(3,i,j)] = vx2l;
         }
         else
         {
            U[c2(0,i,j)] = nr;
            U[c2(1,i,j)] = pr;
            U[c2(2,i,j)] = vx1r;
            U[c2(3,i,j)] = vx2r;
         }
      }
   }
}
