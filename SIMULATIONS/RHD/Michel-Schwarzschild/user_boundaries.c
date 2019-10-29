#include"main.h"

void User_Boundaries(double *B)
{
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
         if(i >= Nx1-gc)
         {
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = velocity_0;
            B(3,i,j) = 0.0;
            B(4,i,j) = 0.0;
         }
      }
   }

#endif
}
