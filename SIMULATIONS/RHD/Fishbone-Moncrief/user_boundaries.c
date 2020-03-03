#include"main.h"

void User_Boundaries(double *B)
{
   double r;
   double t;
   double vr;
   double alpha, beta;
   double M = Black_Hole_Mass;

#if DIM == 1

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for
#endif
   for(int i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
         B(VX1,i) = velocity_0;
      }
   }

#elif DIM == 2

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for collapse(2)
#endif
   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = velocity_0;
            B(3,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for collapse(2)
#endif
   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         r = grid.X1[i];
         t = grid.X2[j];
         alpha = sqrt(r/(r + 2.0*M));
         beta  = 2.0*M/r;
         K_pol = (K - 1.0)*pow(0.01,2.0)/((K - 1.0 - pow(0.01,2.0))*(K));

         if(B(RHO,i,j) < 1.0e-04 || B(VX3,i,j) < 1.0e-05)
         {
            B(RHO,i,j) = 1.0e-04*pow(r/r_in,-1.5);
            B(PRE,i,j) = 1.0e-06*pow(r/r_in,-2.5)*(K - 1.0);
         }

         B(PRE,i,j) = fabs(B(PRE,i,j));
      }
   }

#endif
}
