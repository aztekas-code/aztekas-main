#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_S(double *a, double *uu)
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

   a[0] = 0;
   a[1] = -(((-2*MM*pow(W,2.0)*h*n*pow(x1,2.0))-8*pow(MM,2.0)*pow(W,2.0)*h*n*x1-8*pow(MM,3.0)*pow(W,2.0)*h*n)*pow(w,2.0)+((-2*MM*pow(W,2.0)*h*n*pow(x1,2.0))-8*pow(MM,2.0)*pow(W,2.0)*h*n*x1-8*pow(MM,3.0)*pow(W,2.0)*h*n)*pow(v,2.0)+(2*MM*pow(W,2.0)*h*n*pow(x1,4.0)+2*pow(MM,2.0)*pow(W,2.0)*h*n*pow(x1,3.0))*pow(u,2.0)+MM*pow(W,2.0)*h*n*pow(x1,7.0/2.0)*sqrt(x1+2*MM)*u-2*MM*p*pow(x1,4.0)-10*pow(MM,2.0)*p*pow(x1,3.0)-12*pow(MM,3.0)*p*pow(x1,2.0))/(sqrt(x1)*sqrt(x1+2*MM)*(pow(x1,5.0)+4*MM*pow(x1,4.0)+4*pow(MM,2.0)*pow(x1,3.0)));
   a[2] = -(((-pow(W,2.0)*h*n*pow(x1,2.0))-4*MM*pow(W,2.0)*h*n*x1-4*pow(MM,2.0)*pow(W,2.0)*h*n)*pow(w,2.0)+((-pow(W,2.0)*h*n*pow(x1,2.0))-4*MM*pow(W,2.0)*h*n*x1-4*pow(MM,2.0)*pow(W,2.0)*h*n)*pow(v,2.0)+MM*pow(W,2.0)*h*n*pow(x1,3.0)*pow(u,2.0)+2*MM*pow(W,2.0)*h*n*pow(x1,5.0/2.0)*sqrt(x1+2*MM)*u-2*p*pow(x1,4.0)+(MM*pow(W,2.0)*h*n-8*MM*p)*pow(x1,3.0)+(2*pow(MM,2.0)*pow(W,2.0)*h*n-8*pow(MM,2.0)*p)*pow(x1,2.0))/(pow(x1,5.0)+4*MM*pow(x1,4.0)+4*pow(MM,2.0)*pow(x1,3.0));
   a[3] = 0;
   a[4] = 0;

   return 0;
}
