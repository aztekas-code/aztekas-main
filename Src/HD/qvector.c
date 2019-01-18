#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_Q(double *a, double *uu)
{
   int i;
   double R;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   R = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   for(i = 0; i <= eq; i++)
   {
      if(i == 0)
      {
         a[i] = n;
      }
      else if(i == 1)
      {
         a[i] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
      }
      else if(i == 2)
      {
         a[i] = n*u;
      }
      else if(i == 3)
      {
         a[i] = n*v;
      }
      else if(i == 4)
      {
         a[i] = n*w;
      }
   }

   return 0;
}
