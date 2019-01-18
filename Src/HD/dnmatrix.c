#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
int funct_Dn(double *a, double *uu)
{
   int i;   double R;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   R = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   for(i = 0; i <= 2; i++)
   {
      if(i == 0)
      {
         a[i] = -(sqrt(K*n*p)-n*v)/n;
      }
      else if(i == 1)
      {
         a[i] = (n*v+sqrt(K*n*p))/n;
      }
      else if(i == 2)
      {
         a[i] = v;
      }
   }
     
   return 0;
}
