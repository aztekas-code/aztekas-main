#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_Q2U(double *a, double *uu)
{
   int i, j, k;
   double D, t, m=0, n=0, o=0;
   double R, SS;
   double dWu, dWv, dWw;
   double dhn, dhp;
   double theta, theta_0;
   double f, derf, h, derh, lor;

   R  = sqrt(pow(x2,2.0)+pow(x1,2.0));

   if(dim == 1)
   {
      for(i = 0; i <= Nx1-0; i++)
      {
         D = uu[c1(0,i)];
         t = uu[c1(1,i)];
         m = uu[c1(2,i)];
         n = 0;
         o = 0;

         x1 = X1[i];
         x2 = 0;
         x3 = 0;

         SS = pow(o,2.0)+pow(n,2.0)+pow(m,2.0);

         theta_0 = 0.0;
         f = 2.0;

         while (fabs(f) > 0.0000001)
         {
            h    = 1.0 + (K / (K - 1.0)) * theta_0;
            lor  = sqrt(1.0 + SS / pow(D*h,2.0));
            derh = K / (K - 1.0);

            f    = h * lor - (theta_0 / lor) - (t / D) - 1.0;
            derf = (1.0 / lor) * (derh - 1.0 - theta_0 * ((lor * lor - 1.0) / (lor * lor)) * (derh / h));

            theta = theta_0 - f / derf;
            theta_0 = theta;
         }

         h   = 1.0 + (K / (K - 1.0)) * theta_0;
         lor  = sqrt(1.0 + SS / pow(D*h,2.0));

         a[c1(0,i)] = D / lor;
         a[c1(1,i)] = D*h*lor - t - D;
         a[c1(2,i)] = m / (D * h * lor);
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            D = uu[c2(0,i,j)];
            t = uu[c2(1,i,j)];
            m = uu[c2(2,i,j)];
            n = uu[c2(3,i,j)];
            o = 0;

            x1 = X1[i];
            x2 = X2[j];
            x3 = 0;

            SS = pow(o,2.0)+pow(n,2.0)+pow(m,2.0);

            theta_0 = 0.0;
            f = 2.0;

            while (fabs(f) > 0.000001)
            {
               h    = 1.0 + (K / (K - 1.0)) * theta_0;
               lor  = sqrt(1.0 + SS / pow(D*h,2.0));
               derh = K / (K - 1.0);

               f    = h * lor - (theta_0 / lor) - (t / D) - 1.0;
               derf = (1.0 / lor) * (derh - 1.0 - theta_0 * ((lor * lor - 1.0) / (lor * lor)) * (derh / h));

               theta = theta_0 - f / derf;
               theta_0 = theta;
            }

            h   = 1.0 + (K / (K - 1.0)) * theta;
            lor  = sqrt(1.0 + SS / pow(D*h,2.0));

            a[c2(0,i,j)] = D / lor;
            a[c2(1,i,j)] = D*h*lor - t - D;
            a[c2(2,i,j)] = m / (D * h * lor);
            a[c2(3,i,j)] = n / (D * h * lor);
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
               D = uu[c3(0,i,j,k)];
               t = uu[c3(1,i,j,k)];
               m = uu[c3(2,i,j,k)];
               n = uu[c3(3,i,j,k)];
               o = uu[c3(4,i,j,k)];

               x1 = X1[i];
               x2 = X2[j];
               x3 = X3[k];

               SS = pow(o,2.0)+pow(n,2.0)+pow(m,2.0);

               theta_0 = 0.0;
               f = 2.0;

               while (fabs(f) > 0.000001)
               {
                  h    = 1.0 + (K / (K - 1.0)) * theta_0;
                  lor  = sqrt(1.0 + SS / pow(D*h,2.0));
                  derh = K / (K - 1.0);

                  f    = h * lor - (theta_0 / lor) - (t / D) - 1.0;
                  derf = (1.0 / lor) * (derh - 1.0 - theta_0 * ((lor * lor - 1.0) / (lor * lor)) * (derh / h));

                  theta = theta_0 - f / derf;
                  theta_0 = theta;
               }

               h   = 1.0 + (K / (K - 1.0)) * theta;
               lor  = sqrt(1.0 + SS / pow(D*h,2.0));

               a[c3(0,i,j,k)] = D / lor;
               a[c3(1,i,j,k)] = D*h*lor - t - D;
               a[c3(2,i,j,k)] = m / (D * h * lor);
               a[c3(3,i,j,k)] = n / (D * h * lor);
               a[c3(4,i,j,k)] = o / (D * h * lor);
            }
         }
      }
   }

   return 0;
}
