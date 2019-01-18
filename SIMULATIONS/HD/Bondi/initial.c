/* 
 *  aztekas initial module
 *  Date of creation: 17-01-2019 16:10:43
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

#if dim == 1

   for(i = 0; i <= Nx1; i++)
   {
      U[c1(0,i)] = density_0;
      U[c1(1,i)] = pressure_0;
      U[c1(2,i)] = velocity_0;
   }

#endif
}
