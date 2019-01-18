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

   if (dim == 1)
   {
      x2 = M_PI_2;
   }

   a[0] = pow(x1,2.0)*sin(x2);
   a[1] = 1;
   a[2] = x1;
   a[3] = (x1*sin(x2));

   return 0;
}
