#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_U2Q(double *a, double *uu)
{
   int i, j, k;
   double n, p, u=0, v=0, w=0;
   double R, W, h;
   double dWu, dWv, dWw;
   double dhn, dhp;
   double lapse, dety;

   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         n = uu[c1(0,i)];
         p = uu[c1(1,i)];
         u = uu[c1(2,i)];
         v = 0;
         w = 0;

         x1 = X1[i];
         x2 = 0.0;
         x3 = 0.0;

         R = sqrt(pow(x2,2.0)+pow(x1,2.0));
         W = x1/sqrt(-pow(w,2.0)-pow(x1,2.0)*pow(v,2.0)-pow(x1,2.0)*pow(u,2.0)+pow(x1,2.0));
         h = (K*p+(K-1)*n)/((K-1)*n);

         lapse = 1;
         dety  = x1;

         a[c1(0,i)] = n*W;
         a[c1(1,i)] = h*n*pow(W,2.0)-n*W-p;
         a[c1(2,i)] = h*n*u*pow(W,2.0);
         a[c1(3,i)] = h*n*v*pow(W,2.0);
         a[c1(4,i)] = h*n*w*pow(W,2.0);
      }
   }
   else if(dim == 2)
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

            x1 = X1[i];
            x2 = X2[j];
            x3 = 0.0;

            R = sqrt(pow(x2,2.0)+pow(x1,2.0));
            W = x1/sqrt(-pow(w,2.0)-pow(x1,2.0)*pow(v,2.0)-pow(x1,2.0)*pow(u,2.0)+pow(x1,2.0));
            h = (K*p+(K-1)*n)/((K-1)*n);

            lapse = 1;
            dety  = x1;

            a[c2(0,i,j)] = n*W;
            a[c2(1,i,j)] = h*n*pow(W,2.0)-n*W-p;
            a[c2(2,i,j)] = h*n*u*pow(W,2.0);
            a[c2(3,i,j)] = h*n*v*pow(W,2.0);
            a[c2(4,i,j)] = h*n*w*pow(W,2.0);
         }
      }
   }
   if(dim == 3)
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

               x1 = X1[i];
               x2 = X2[j];
               x3 = X3[k];

               R = sqrt(pow(x2,2.0)+pow(x1,2.0));
               W = x1/sqrt(-pow(w,2.0)-pow(x1,2.0)*pow(v,2.0)-pow(x1,2.0)*pow(u,2.0)+pow(x1,2.0));
               h = (K*p+(K-1)*n)/((K-1)*n);

               lapse = 1;
               dety  = 1;

               a[c3(0,i,j,k)] = n*W;
               a[c3(1,i,j,k)] = h*n*pow(W,2.0)-n*W-p;
               a[c3(2,i,j,k)] = h*n*u*pow(W,2.0);
               a[c3(3,i,j,k)] = h*n*v*pow(W,2.0);
               a[c3(4,i,j,k)] = h*n*w*pow(W,2.0);
            }
         }
      }
   }

   return 0;
}
