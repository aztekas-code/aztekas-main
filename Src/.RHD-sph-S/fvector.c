#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_F(double *a, double *uu)
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
   W = (x1*(sin(x2)))/sqrt((-pow(w,2.0))-pow(sin(x2),2.0)*pow(v,2.0)+(2*MM*x1-pow(x1,2.0))*pow(sin(x2),2.0)*pow(u,2.0)+pow(x1,2.0)*pow(sin(x2),2.0));
   h = (K*p+(K-1)*n)/((K-1)*n);

   a[0] = ((W*n*x1-2*MM*W*n)*u)/x1;
   a[1] = (((pow(W,2.0)*h-W)*n*x1+(2*MM*W-2*MM*pow(W,2.0)*h)*n)*u)/x1;
   a[2] = ((pow(W,2.0)*h*n*x1-2*MM*pow(W,2.0)*h*n)*pow(u,2.0)+p*x1)/x1;
   a[3] = ((pow(W,2.0)*h*n*x1-2*MM*pow(W,2.0)*h*n)*u*v)/x1;
   a[4] = ((pow(W,2.0)*h*n*x1-2*MM*pow(W,2.0)*h*n)*u*w)/x1;

   return 0;
}
