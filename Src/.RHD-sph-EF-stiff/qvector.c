#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_Q(double *a, double *uu)
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

   R = x1;
   W = (x1*(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
   h = (100000*K*p+(K-1)*n)/((100000*K-100000)*n);

   a[0] = W*n;
   a[1] = (pow(W,2.0)*h-W)*n-p;
   a[2] = pow(W,2.0)*h*n*u;
   a[3] = pow(W,2.0)*h*n*v;
   a[4] = pow(W,2.0)*h*n*w;

   return 0;
}
