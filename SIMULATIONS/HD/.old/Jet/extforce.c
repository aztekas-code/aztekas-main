/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 18:32:00
 *  author: Alejandro Aguayo Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

int EXTFORCE(double *a, double *uu)
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

   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = 0;
      }
      else if(i == 1)
      {
         a[i] = 0;
      }
      else if(i == 2)
      {
         a[i] = 0;
      }
      else if(i == 3)
      {
         a[i] = 0;
      }
      else if(i == 4)
      {
         a[i] = 0;
      }
   }

   return 0;
}
