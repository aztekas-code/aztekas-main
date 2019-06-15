#include"main.h"
    
void Prim2Cons_All(double *a, double *uu)
{
   int i, j, k;
   double R;
   double n, p, u, v, w;
   R = sqrt(x1*x1 + x2*x2);
   
#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      n = uu(0,i);
      p = uu(1,i);
      u = uu(2,i);
      v = 0;
      w = 0;
 
      a(0,i) = n;
      a(1,i) = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
      a(2,i) = n*u;
      a(3,i) = n*v;
      a(4,i) = n*w;
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         n = uu(0,i,j);
         p = uu(1,i,j);
         u = uu(2,i,j);
         v = uu(3,i,j);
         w = 0;
 
         a(0,i,j) = n;
         a(1,i,j) = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
         a(2,i,j) = n*u;
         a(3,i,j) = n*v;
         a(4,i,j) = n*w;
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         n = uu(0,i,j);
         p = uu(1,i,j);
         u = uu(2,i,j);
         v = uu(3,i,j);
         w = uu(4,i,j);
 
         a(0,i,j) = n;
         a(1,i,j) = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
         a(2,i,j) = n*u;
         a(3,i,j) = n*v;
         a(4,i,j) = n*w;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            n = uu(0,i,j,k);
            p = uu(1,i,j,k);
            u = uu(2,i,j,k);
            v = uu(3,i,j,k);
            w = uu(4,i,j,k);
 
            a(0,i,j,k) = n;
            a(1,i,j,k) = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
            a(2,i,j,k) = n*u;
            a(3,i,j,k) = n*v;
            a(4,i,j,k) = n*w;
         }
      }
   }

#endif
}
