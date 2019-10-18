#include"main.h"
    
int Cons2Prim(double *u, double *q)
{
   int i, j, k;
   double D, E, S1, S2, S3;
   
#if DIM == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      D  = q(RHO,i);
      E  = q(PRE,i);
      S1 = q(VX1,i);
      S2 = 0;
      S3 = 0;
 
      u(RHO,i) = D;
      #if EOS == IDEAL
      u(PRE,i) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
      #elif EOS == DUST
      u(PRE,i) = 0.0;
      #endif
      u(VX1,i) = S1/D;
   }

#elif DIM == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = q(RHO,i,j);
         E  = q(PRE,i,j);
         S1 = q(VX1,i,j);
         S2 = q(VX2,i,j);
         S3 = 0;
 
         u(RHO,i,j) = D;
         #if EOS == IDEAL
         u(PRE,i,j) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #elif EOS == DUST
         u(PRE,i,j) = 0.0;
         #endif
         u(VX1,i,j) = S1/D;
         u(VX2,i,j) = S2/D;
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         D  = q(RHO,i,j);
         E  = q(PRE,i,j);
         S1 = q(VX1,i,j);
         S2 = q(VX2,i,j);
         S3 = q(VX3,i,j);
 
         u(RHO,i,j) = D;
         #if EOS == IDEAL
         u(PRE,i,j) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #elif EOS == DUST
         u(PRE,i,j) = 0.0;
         #endif
         u(VX1,i,j) = S1/D;
         u(VX2,i,j) = S2/D;
         u(VX3,i,j) = S3/D;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            D  = q(0,i,j,k);
            E  = q(1,i,j,k);
            S1 = q(2,i,j,k);
            S2 = q(3,i,j,k);
            S3 = q(4,i,j,k);
 
            u(0,i,j,k) = D;
            #if EOS == IDEAL
            u(1,i,j,k) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
            #elif EOS == DUST
            u(1,i,j,k) = 0.0;
            #endif
            u(2,i,j,k) = S1/D;
            u(3,i,j,k) = S2/D;
            u(4,i,j,k) = S3/D;
         }
      }
   }

#endif
   
   return 0;
}
