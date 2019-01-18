#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_S(double *a, double *uu)
{
   int i;
   double r;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim == 1){x2 = M_PI_2;};
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
   double E = 0.5*n*(u*u + v*v + w*w) + p/(K-1);
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = -2.0*u*n/x1 - v*n*cos(x2)/(x1*sin(x2));
      }
      else if(i == 1)
      {
         a[i] = -2.0*u*(E+p)/x1 - v*(E+p)*cos(x2)/(x1*sin(x2)) -(n*u)/pow(x1,2.0);
      }
      else if(i == 2)
      {
         a[i] = -2.0*u*(n*u)/x1 - v*(n*u)*cos(x2)/(x1*sin(x2)) -n/pow(x1,2.0);
      }
      else if(i == 3)
      {
         a[i] = -2.0*u*(n*v)/x1 - v*(n*v)*cos(x2)/(x1*sin(x2));
      }
      else if(i == 4)
      {
         a[i] = 0;
      }
   }

   return 0;
}
