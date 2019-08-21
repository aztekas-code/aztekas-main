/* 
 *  aztekas initial module
 *  Date of creation: 17-05-2019 19:07:17
 *  author: Alejandro Aguayo-Ortiz 
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

#if dim == 1
   for(i = 0; i <= Nx1; i++)
   {
      U[c1(0,i)] = density_0;
      U[c1(1,i)] = pressure_0;
      U[c1(2,i)] = velocity_0;
   }

#elif dim == 2

   for(i = 0; i <= Nx1; i++)
   {   
      for(j = 0; j <= Nx2; j++)
      {
         U[c2(0,i,j)] = density_0;
         U[c2(1,i,j)] = pressure_0;
         U[c2(2,i,j)] = velocity_0;
         U[c2(3,i,j)] = 0.0;
      }
   }

#endif
}
