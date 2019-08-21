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

   a[0] = x1;
   a[1] = 1.0;
   a[2] = 1.0;
   a[3] = x1;

   return 0;
}
