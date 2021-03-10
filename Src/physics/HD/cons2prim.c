/**
 * @file /HD/cons2prim.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Computes the value of the conservative variables \f$ Q \f$
 * from the primitive variables in HD.
 */

#include"main.h"
    
int Cons2Prim(double *u, double *q)
{
   double D, E, S1, S2, S3;
   double CONS[3];
   eos_ eos;
   gauge_ local_grid;
   
#if DIM == 1

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,E,S1,S2,S3,CONS,eos,local_grid)
#endif
   for(int i = gc; i <= Nx1-gc; i++)
   {
      D  = q(0,i);
      E  = q(1,i);
      S1 = q(2,i);
      S2 = 0;
      S3 = 0;

      u(RHO,i) = D;
      #if EOS == IDEAL
      u(PRE,i) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
      #elif EOS == DUST
      u(PRE,i) = 0.0;
      #elif EOS == HELMHOLTZ
      CONS[0] = D;
      CONS[1] = 0.0;
      CONS[2] = (E - 0.5 * (S1 * S1 + S2 * S2 + S3 * S3) / D);

      EoS(&eos,CONS,&local_grid);

      u(PRE,i) = eos.p;
      #endif
      u(VX1,i) = S1/D;
   }

#elif DIM == 2

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,E,S1,S2,S3) collapse(2)
#endif   
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         D  = q(RHO,i,j);
         E  = q(ENE,i,j);
         S1 = q(MX1,i,j);
         S2 = q(MX2,i,j);
         S3 = 0;
 
         u(RHO,i,j) = D;
         #if EOS == IDEAL
         u(PRE,i,j) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #elif EOS == DUST
         u(PRE,i,j) = 0.0;
         #elif EOS == HELMHOLTZ
         CONS[0] = D;
         CONS[1] = 0.0;
         CONS[2] = E - 0.5 * (S1 * S1 + S2 * S2 + S3 * S3) / D;

         EoS(&eos,CONS,&local_grid);

         u(PRE,i,j) = eos.p;
         #endif
         u(VX1,i,j) = S1/D;
         u(VX2,i,j) = S2/D;
      }
   }

#elif DIM == 4

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,E,S1,S2,S3) collapse(2)
#endif   
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         D  = q(RHO,i,j);
         E  = q(ENE,i,j);
         S1 = q(MX1,i,j);
         S2 = q(MX2,i,j);
         S3 = q(MX3,i,j);
 
         u(RHO,i,j) = D;
         #if EOS == IDEAL
         u(PRE,i,j) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
         #elif EOS == DUST
         u(PRE,i,j) = 0.0;
         #elif EOS == HELMHOLTZ
         CONS[0] = D;
         CONS[1] = 0.0;
         CONS[2] = E - 0.5 * (S1 * S1 + S2 * S2 + S3 * S3) / D;

         EoS(&eos,CONS,&local_grid);

         u(PRE,i,j) = eos.p;
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
            D  = q(RHO,i,j,k);
            E  = q(ENE,i,j,k);
            S1 = q(MX1,i,j,k);
            S2 = q(MX2,i,j,k);
            S3 = q(MX3,i,j,k);
 
            u(RHO,i,j,k) = D;
            #if EOS == IDEAL
            u(PRE,i,j,k) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
            #elif EOS == DUST
            u(PRE,i,j,k) = 0.0;
            #elif EOS == HELMHOLTZ
            CONS[0] = D;
            CONS[1] = 0.0;
            CONS[2] = E - 0.5 * (S1 * S1 + S2 * S2 + S3 * S3) / D;

            EoS(&eos,CONS,&local_grid);

            u(PRE,i,j) = eos.p;
            #endif
            u(VX1,i,j,k) = S1/D;
            u(VX2,i,j,k) = S2/D;
            u(VX3,i,j,k) = S3/D;
         }
      }
   }

#endif
   
   return 0;
}
