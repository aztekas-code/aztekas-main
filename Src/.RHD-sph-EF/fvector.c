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
   W = x1/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(v,2.0)+pow(x1,3.0)*pow(u,2.0)-pow(x1,3.0)-2*MM*pow(x1,2.0))/(x1+2*MM));
   h = (K*p+(K-1)*n)/((K-1)*n);

   a[0] = (W*n*pow(x1,3.0/2.0)*sqrt(x1+2*MM)*u-2*MM*W*n*x1-4*pow(MM,2.0)*W*n)/(sqrt(x1)*pow(2*MM+x1,1.5));
   a[1] = ((pow(W,2.0)*h-W)*n*pow(x1,3.0/2.0)*sqrt(x1+2*MM)*u+(2*MM*p+(2*MM*W-2*MM*pow(W,2.0)*h)*n)*x1+4*pow(MM,2.0)*p+(4*pow(MM,2.0)*W-4*pow(MM,2.0)*pow(W,2.0)*h)*n)/(sqrt(x1)*pow(2*MM+x1,1.5));
   a[2] = (sqrt(x1)*sqrt(x1+2*MM)*(pow(W,2.0)*h*n*x1*pow(u,2.0)+p*x1+2*MM*p)+((-2*MM*pow(W,2.0)*h*n*x1)-4*pow(MM,2.0)*pow(W,2.0)*h*n)*u)/(sqrt(x1)*pow(2*MM+x1,1.5));
   a[3] = (pow(W,2.0)*h*n*pow(x1,3.0/2.0)*sqrt(x1+2*MM)*u*v+((-2*MM*pow(W,2.0)*h*n*x1)-4*pow(MM,2.0)*pow(W,2.0)*h*n)*v)/(sqrt(x1)*pow(2*MM+x1,1.5));
   a[4] = (pow(W,2.0)*h*n*pow(x1,3.0/2.0)*sqrt(x1+2*MM)*u*w+((-2*MM*pow(W,2.0)*h*n*x1)-4*pow(MM,2.0)*pow(W,2.0)*h*n)*w)/(sqrt(x1)*pow(2*MM+x1,1.5));

   return 0;
}
