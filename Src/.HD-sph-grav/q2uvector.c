#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"
    
int funct_Q2U(double *a, double *uu)
{
   int i, j, k;
   double r;
   double D, E, m1, m2, m3;
   r = sqrt(x1*x1 + x2*x2);
   
   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         D  = uu[c1(0,i)];
         E  = uu[c1(1,i)];
         m1 = uu[c1(2,i)];
         m2 = 0;
         m3 = 0;
    
         a[c1(0,i)] = D;
         a[c1(1,i)] = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
         a[c1(2,i)] = m1/D;
         a[c1(3,i)] = m2/D;
         a[c1(4,i)] = m3/D;
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            D  = uu[c2(0,i,j)];
            E  = uu[c2(1,i,j)];
            m1 = uu[c2(2,i,j)];
            m2 = uu[c2(3,i,j)];
            m3 = 0;
    
            a[c2(0,i,j)] = D;
            a[c2(1,i,j)] = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
            a[c2(2,i,j)] = m1/D;
            a[c2(3,i,j)] = m2/D;
            a[c2(4,i,j)] = m3/D;
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
               D  = uu[c3(0,i,j,k)];
               E  = uu[c3(1,i,j,k)];
               m1 = uu[c3(2,i,j,k)];
               m2 = uu[c3(3,i,j,k)];
               m3 = uu[c3(4,i,j,k)];
    
               a[c3(0,i,j,k)] = D;
               a[c3(1,i,j,k)] = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
               a[c3(2,i,j,k)] = m1/D;
               a[c3(3,i,j,k)] = m2/D;
               a[c3(4,i,j,k)] = m3/D;
            }
         }
      }
   }
   
   return 0;
}
