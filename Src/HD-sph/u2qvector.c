#include"main.h"
    
void Prim2Cons_All(double *a, double *uu)
{
   int i, j, k;
   double R;
   double n, p, u, v, w;
   R = sqrt(x1*x1 + x2*x2);
   
   if(DIM == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         n = uu[c1(0,i)];
         p = uu[c1(1,i)];
         u = uu[c1(2,i)];
         v = 0;
         w = 0;
    
         a[c1(0,i)] = n;
         a[c1(1,i)] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
         a[c1(2,i)] = n*u;
         a[c1(3,i)] = n*v;
         a[c1(4,i)] = n*w;
      }
   }
   else if(DIM == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            n = uu[c2(0,i,j)];
            p = uu[c2(1,i,j)];
            u = uu[c2(2,i,j)];
            v = uu[c2(3,i,j)];
            w = 0;
    
            a[c2(0,i,j)] = n;
            a[c2(1,i,j)] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
            a[c2(2,i,j)] = n*u;
            a[c2(3,i,j)] = n*v;
            a[c2(4,i,j)] = n*w;
         }
      }
   }
   else if(DIM == 4)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            n = uu[c2(0,i,j)];
            p = uu[c2(1,i,j)];
            u = uu[c2(2,i,j)];
            v = uu[c2(3,i,j)];
            w = uu[c2(4,i,j)];
    
            a[c2(0,i,j)] = n;
            a[c2(1,i,j)] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
            a[c2(2,i,j)] = n*u;
            a[c2(3,i,j)] = n*v;
            a[c2(4,i,j)] = n*w;
         }
      }
   }
   if(DIM == 3)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            for(k = 0; k <= Nx3; k++)
            {
               n = uu[c3(0,i,j,k)];
               p = uu[c3(1,i,j,k)];
               u = uu[c3(2,i,j,k)];
               v = uu[c3(3,i,j,k)];
               w = uu[c3(4,i,j,k)];
    
               a[c3(0,i,j,k)] = n;
               a[c3(1,i,j,k)] = ((K-1)*n*pow(w,2.0)+(K-1)*n*pow(v,2.0)+(K-1)*n*pow(u,2.0)+2*p)/(2*K-2);
               a[c3(2,i,j,k)] = n*u;
               a[c3(3,i,j,k)] = n*v;
               a[c3(4,i,j,k)] = n*w;
            }
         }
      }
   }
}
