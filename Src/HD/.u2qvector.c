#include"../Headers/main.h"
    
int funct_U2Q(double *a, double *u)
{
   int n, i, j, k;
   double E, e;
   double p[eq+1];
   double x[4];
   
#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      x[0] = time;
      x[1] = X1[i];
      x[2] = 0.0;
      x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      x[2] = M_PI_2;
      #endif

      for(n = 0; n < eq; n++)
      {
         p[n] = u[c1(n,i)];
      }

      #if EOS == IDEAL  
      EoS_Ideal(&e,p,x);
      #endif
      
      a[c1(0,i)] = p[0];
      a[c1(1,i)] = 0.5 * p[0]*(p[2]*p[2] + p[3]*p[3] + p[4]*p[4]) + p[0]*e;
      a[c1(2,i)] = p[0]*p[2];
      a[c1(3,i)] = p[0]*p[3];
      a[c1(4,i)] = p[0]*p[4];
   }

#elif DIM == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         x[0] = time;
         x[1] = X1[i];
         x[2] = 0.0;
         x[3] = 0.0;
         #if POLAR == TRUE
         x[2] = M_PI_2;
         #endif
       
         for(n = 0; n < eq; n++)
         {
            p[n] = u[c2(n,i,j)];
         }
       
         #if EOS == IDEAL  
         EoS_Ideal(&e,p,x);
         #endif
      
         a[c2(0,i,j)] = p[0];
         a[c2(1,i,j)] = 0.5 * p[0]*(p[2]*p[2] + p[3]*p[3] + p[4]*p[4]) + p[0]*e;
         a[c2(2,i,j)] = p[0]*p[2];
         a[c2(3,i,j)] = p[0]*p[3];
         a[c2(4,i,j)] = p[0]*p[4];
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         n = u[c2(0,i,j)];
         p = u[c2(1,i,j)];
         u = u[c2(2,i,j)];
         v = u[c2(3,i,j)];
         w = u[c2(4,i,j)];
 
         a[c2(0,i,j)] = n;
         a[c2(1,i,j)] = ((K-1.0)*n*pow(w,2.0)+(K-1.0)*n*pow(v,2.0)+(K-1.0)*n*pow(u,2.0)+2.0*p)/(2.0*K-2.0);
         a[c2(2,i,j)] = n*u;
         a[c2(3,i,j)] = n*v;
         a[c2(4,i,j)] = n*w;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            n = u[c3(0,i,j,k)];
            p = u[c3(1,i,j,k)];
            u = u[c3(2,i,j,k)];
            v = u[c3(3,i,j,k)];
            w = u[c3(4,i,j,k)];
 
            a[c3(0,i,j,k)] = n;
            a[c3(1,i,j,k)] = ((K-1.0)*n*pow(w,2.0)+(K-1.0)*n*pow(v,2.0)+(K-1.0)*n*pow(u,2.0)+2.0*p)/(2.0*K-2.0);
            a[c3(2,i,j,k)] = n*u;
            a[c3(3,i,j,k)] = n*v;
            a[c3(4,i,j,k)] = n*w;
         }
      }
   }

#endif

   return 0;
}
