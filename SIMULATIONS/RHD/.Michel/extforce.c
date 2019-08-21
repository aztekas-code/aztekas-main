/* 
 *  aztekas boundaries module
 *  Date of creation: 17-05-2019 19:07:17
 *  author: Alejandro Aguayo-Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

void EXTFORCE(double *a, double *uu)
{
   int i;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
#if dim == 1
   u = uu[2];
#elif dim == 2
   u = uu[2];
   v = uu[3];
#elif dim == 3 || dim == 4
   u = uu[2];
   v = uu[3];
   w = uu[4];
#endif

   a[0] = 0;
   a[1] = 0;
   a[2] = 0;
   a[3] = 0;
   a[4] = 0;
}
