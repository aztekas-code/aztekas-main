#include"main.h"
    
int funct_Q2U(double *a, double *uu)
{
   int i, j, k;
   double D, E, m1, m2, m3;
   
#if dim == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      D  = uu[c1(0,i)];
      E  = uu[c1(1,i)];
      m1 = uu[c1(2,i)];
      m2 = 0;
      m3 = 0;
 
      a[c1(0,i)] = D;
      a[c1(1,i)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(m3,2.0)+(1.0-K)*pow(m2,2.0)+(1.0-K)*pow(m1,2.0))/(2.0*D);
      a[c1(2,i)] = m1/D;
      a[c1(3,i)] = m2/D;
      a[c1(4,i)] = m3/D;
   }

#elif dim == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = uu[c2(0,i,j)];
         E  = uu[c2(1,i,j)];
         m1 = uu[c2(2,i,j)];
         m2 = uu[c2(3,i,j)];
         m3 = 0;
 
         a[c2(0,i,j)] = D;
         a[c2(1,i,j)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(m3,2.0)+(1.0-K)*pow(m2,2.0)+(1.0-K)*pow(m1,2.0))/(2.0*D);
         a[c2(2,i,j)] = m1/D;
         a[c2(3,i,j)] = m2/D;
         a[c2(4,i,j)] = m3/D;
      }
   }

#elif dim == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = uu[c2(0,i,j)];
         E  = uu[c2(1,i,j)];
         m1 = uu[c2(2,i,j)];
         m2 = uu[c2(3,i,j)];
         m3 = uu[c2(4,i,j)];
 
         a[c2(0,i,j)] = D;
         a[c2(1,i,j)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(m3,2.0)+(1.0-K)*pow(m2,2.0)+(1.0-K)*pow(m1,2.0))/(2.0*D);
         a[c2(2,i,j)] = m1/D;
         a[c2(3,i,j)] = m2/D;
         a[c2(4,i,j)] = m3/D;
      }
   }

#elif dim == 3

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
            a[c3(1,i,j,k)] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(m3,2.0)+(1.0-K)*pow(m2,2.0)+(1.0-K)*pow(m1,2.0))/(2.0*D);
            a[c3(2,i,j,k)] = m1/D;
            a[c3(3,i,j,k)] = m2/D;
            a[c3(4,i,j,k)] = m3/D;
         }
      }
   }

#endif
   
   return 0;
}
