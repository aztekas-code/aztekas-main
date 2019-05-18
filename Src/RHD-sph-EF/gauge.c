#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int GAUGE(double *a, double g1, double g2, double g3)
{
   x1 = g1;
   x2 = g2;
   x3 = g3;

   a[0] = sqrt(x1)/sqrt(x1+2*MM);
   a[1] = pow(x1,3.0/2.0)*sqrt(x1+2*MM);
   a[2] = pow(x1,2.0);

   return 0;
}
