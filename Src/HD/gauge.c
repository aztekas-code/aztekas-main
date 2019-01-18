#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int GAUGE(double *a, double g1, double g2, double g3)
{
   x1 = g1;
   x2 = g2;
   x3 = g3;

   if(alfa == 0)
   {
      a[0] = 1;
      a[1] = 1;
      a[2] = 1;
   }
   else if(alfa == 1)
   {
      a[0] = 1;
      a[1] = x1;
      a[2] = x1;
   }
   else if(alfa == 2)
   {
      if(dim == 1)
      {
         x2 = M_PI_2;
      }

      a[0] = 1;
      a[1] = pow(x1*sin(x2),2.0);
      a[2] = pow(x1*sin(x2),2.0);
   }
   else
   {
      printf("alfa parameter must be 0, 1 or 2");
      exit(1);
   }

   return 0;
}
