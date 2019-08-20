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
   W = (x1*(sin(x2)))/sqrt((-pow(w,2.0))-pow(sin(x2),2.0)*pow(v,2.0)+(2*MM*x1-pow(x1,2.0))*pow(sin(x2),2.0)*pow(u,2.0)+pow(x1,2.0)*pow(sin(x2),2.0));
   h = (K*p+(K-1)*n)/((K-1)*n);

   a[0] = 0;
   a[1] = -(MM*pow(W,2.0)*h*n*u)/pow(x1,2.0);
   a[2] = ((pow(W,2.0)*h*n*x1-2*MM*pow(W,2.0)*h*n)*pow(w,2.0)+(pow(W,2.0)*h*n*x1-2*MM*pow(W,2.0)*h*n)*pow(sin(x2),2.0)*pow(v,2.0)+(2*pow(MM,2.0)*pow(W,2.0)*h*n*x1-MM*pow(W,2.0)*h*n*pow(x1,2.0))*pow(sin(x2),2.0)*pow(u,2.0)+(2*p*pow(x1,3.0)+((-4*MM*p)-MM*pow(W,2.0)*h*n)*pow(x1,2.0))*pow(sin(x2),2.0))/((pow(x1,4.0)-2*MM*pow(x1,3.0))*pow(sin(x2),2.0));
   a[3] = (pow(W,2.0)*h*n*cos(x2)*pow(w,2.0)+p*pow(x1,2.0)*cos(x2)*pow(sin(x2),2.0))/(pow(x1,2.0)*pow(sin(x2),3.0));
   a[4] = 0;

   return 0;
}
