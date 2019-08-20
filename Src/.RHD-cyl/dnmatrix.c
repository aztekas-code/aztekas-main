#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_Dn(double *a, double *uu)
{
   int i;
   double n, p, u=0, v=0, w=0;
   double R, V;
   double yt22, Bt2, Vt2, c;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}

   R    = sqrt(pow(x2,2.0)+pow(x1,2.0));
   V    = (pow(w,2.0)+pow(x1,2.0)*pow(v,2.0)+pow(x1,2.0)*pow(u,2.0))/pow(x1,2.0);
   yt22 = 1;
   Bt2  = 0;
   Vt2  = v;
   c    = sqrt((pow(K,2.0)-K)*p/(K*p+(K-1)*n));

   a[0] = -((c)*sqrt(pow(c,2.0)*yt22*pow(V,2.0)+((-pow(c,2.0)-1)*yt22+(1-pow(c,2.0))*pow(Vt2,2.0))*V+yt22+(pow(c,2.0)-1)*pow(Vt2,2.0))+Bt2*pow(c,2.0)*V+(1-pow(c,2.0))*Vt2-Bt2)/(pow(c,2.0)*V-1);
   a[1] = ((c)*sqrt(pow(c,2.0)*yt22*pow(V,2.0)+((-pow(c,2.0)-1)*yt22+(1-pow(c,2.0))*pow(Vt2,2.0))*V+yt22+(pow(c,2.0)-1)*pow(Vt2,2.0))-Bt2*pow(c,2.0)*V+(pow(c,2.0)-1)*Vt2+Bt2)/(pow(c,2.0)*V-1);
   a[2] = Vt2-Bt2;

   return 0;
}
