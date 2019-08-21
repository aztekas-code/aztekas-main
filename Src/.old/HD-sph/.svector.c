#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
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

   double E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);
   
   a[0] = -(2*n*sin(x2)*u*x1+n*x1*cos(x2)*v)/(x1*sin(x2)*x1);
   a[1] = -((2*sin(x2)*u*x1+x1*cos(x2)*v)*E+2*p*sin(x2)*u*x1+p*x1*cos(x2)*v)/(x1*sin(x2)*x1);
   a[2] = ((n*sin(x2)*pow(w,2.0)+n*sin(x2)*pow(v,2.0)-2*n*sin(x2)*pow(u,2.0))*x1-n*x1*cos(x2)*u*v)/(x1*sin(x2)*x1);
   a[3] = ((n*cos(x2)*pow(w,2.0)-3*n*sin(x2)*u*v)*x1-n*x1*cos(x2)*pow(v,2.0))/(x1*sin(x2)*x1);
   a[4] = -((n*cos(x2)*v+3*n*sin(x2)*u)*w*x1+n*x1*cos(x2)*v*w)/(x1*sin(x2)*x1);

   return 0;
}
