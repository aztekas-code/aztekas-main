#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   double E;

   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim == 1){x2 = M_PI_2;};
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
   E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = (-(n*cos(x2)*v)/(sin(x2)*x1))-(2.0*n*u)/x1;
      }
      else if(i == 1)
      {
         a[i] = (-((p+E)*cos(x2)*v)/(sin(x2)*x1))-(2.0*(p+E)*u)/x1-(n*u)/pow(x1,2.0);
      }
      else if(i == 2)
      {
         a[i] = (-(n*cos(x2)*u*v)/(sin(x2)*x1))+(n*(pow(w,2.0)+pow(v,2.0)))/x1-(2.0*n*pow(u,2.0))/x1-n/pow(x1,2.0);
      }
      else if(i == 3)
      {
         a[i] = (-(n*cos(x2)*pow(v,2.0))/(sin(x2)*x1))+(n*cos(x2)*pow(w,2.0))/(x1*sin(x2))-(3.0*n*u*v)/x1;
      }
      else if(i == 4)
      {
         a[i] = (-(n*cos(x2)*v*w)/(sin(x2)*x1))-(n*cos(x2)*v*w)/(x1*sin(x2))-(3.0*n*u*w)/x1;
      }
   }

   return 0;
}
