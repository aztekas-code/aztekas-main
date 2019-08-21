#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_G(double *a, double *uu)
{
   int i;
   double n, p, u=0, v=0, w=0;
   double R, W, h;
   double dWu, dWv, dWw;
   double dhn, dhp;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}

   R = sqrt(pow(x2,2.0)+pow(x1,2.0));
   W = 1/sqrt((-pow(w,2.0))-pow(v,2.0)-pow(u,2.0)+1);
   h = (K*p+(K-1)*n)/((K-1)*n);

   a[0] = W*n*v;
   a[1] = (pow(W,2.0)*h-W)*n*v;
   a[2] = pow(W,2.0)*h*n*u*v;
   a[3] = pow(W,2.0)*h*n*pow(v,2.0)+p;
   a[4] = pow(W,2.0)*h*n*v*w;

   return 0;
}
