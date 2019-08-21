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

         R = x1;
         W = (x1*(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
         h = (K*p+(K-1)*n)/((K-1)*n);

         lapse = sqrt(x1)/sqrt(x1+2*MM);
         dety  = pow(x1,3.0/2.0)*sqrt(x1+2*MM)*(sin(x2));

         a[c1(0,i)] = W*n;
         a[c1(1,i)] = (pow(W,2.0)*h-W)*n-p;
         a[c1(2,i)] = pow(W,2.0)*h*n*u;
         a[c1(3,i)] = pow(W,2.0)*h*n*v;
         a[c1(4,i)] = pow(W,2.0)*h*n*w;
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

            R = x1;
            W = (x1*(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
            h = (K*p+(K-1)*n)/((K-1)*n);

            lapse = sqrt(x1)/sqrt(x1+2*MM);
            dety  = pow(x1,3.0/2.0)*sqrt(x1+2*MM)*(sin(x2));

            a[c2(0,i,j)] = W*n;
            a[c2(1,i,j)] = (pow(W,2.0)*h-W)*n-p;
            a[c2(2,i,j)] = pow(W,2.0)*h*n*u;
            a[c2(3,i,j)] = pow(W,2.0)*h*n*v;
            a[c2(4,i,j)] = pow(W,2.0)*h*n*w;
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

               R = x1;
               W = (x1*(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
               h = (K*p+(K-1)*n)/((K-1)*n);

               lapse = sqrt(x1)/sqrt(x1+2*MM);
               dety  = sqrt(x1)/sqrt(x1+2*MM);

               a[c3(0,i,j,k)] = W*n;
               a[c3(1,i,j,k)] = (pow(W,2.0)*h-W)*n-p;
               a[c3(2,i,j,k)] = pow(W,2.0)*h*n*u;
               a[c3(3,i,j,k)] = pow(W,2.0)*h*n*v;
               a[c3(4,i,j,k)] = pow(W,2.0)*h*n*w;
            }
         }
      }
   }

   return 0;
}
