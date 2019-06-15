#include"main.h"
    
int funct_Q2U(double *a, double *uu)
{
   int i, j, k;
   double R;
   double D, E, m1, m2, m3;
   R = sqrt(x1*x1 + x2*x2);
   
#if DIM == 1
{
   for(i = 0; i <= Nx1; i++)
   {
      D  = uu(0,i);
      E  = uu(1,i);
      m1 = uu(2,i);
      m2 = 0;
      m3 = 0;
 
      a(0,i) = D;
      a(1,i) = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
      a(2,i) = m1/D;
      a(3,i) = m2/D;
      a(4,i) = m3/D;
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         D  = uu(0,i,j);
         E  = uu(1,i,j);
         m1 = uu(2,i,j);
         m2 = uu(3,i,j);
         m3 = 0;
 
         a(0,i,j) = D;
         a(1,i,j) = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
         a(2,i,j) = m1/D;
         a(3,i,j) = m2/D;
         a(4,i,j) = m3/D;
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         D  = uu(0,i,j);
         E  = uu(1,i,j);
         m1 = uu(2,i,j);
         m2 = uu(3,i,j);
         m3 = uu(4,i,j);
 
         a(0,i,j) = D;
         a(1,i,j) = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
         a(2,i,j) = m1/D;
         a(3,i,j) = m2/D;
         a(4,i,j) = m3/D;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            D  = uu(0,i,j,k);
            E  = uu(1,i,j,k);
            m1 = uu(2,i,j,k);
            m2 = uu(3,i,j,k);
            m3 = uu(4,i,j,k);
 
            a(0,i,j,k) = D;
            a(1,i,j,k) = ((2*K-2)*D*E+(1-K)*pow(m3,2.0)+(1-K)*pow(m2,2.0)+(1-K)*pow(m1,2.0))/(2*D);
            a(2,i,j,k) = m1/D;
            a(3,i,j,k) = m2/D;
            a(4,i,j,k) = m3/D;
         }
      }
   }

#endif
   
   return 0;
}
