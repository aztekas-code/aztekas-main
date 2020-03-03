#include"main.h"

void User_Boundaries(double *B)
{
   FILE *file;
   int idum;
   double dum, dumrho, dumpre;
   char line[100];
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
         //B(VX1,i) = velocity_0;
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

   for(int j = 0; j <= Nx2; j++)
   {
      file = fopen("./Michel/analytic","r");
      for(int i = 0; i <= Nx1; i++)
      {
         r = grid.X1[i];
         t = grid.X2[j];
         alpha = sqrt(r/(r + 2.0*M));
         beta  = 2.0*M/r;
         vr = (r/(r + 2.0*M))*B(VX1,i,j);

         idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,&dumrho,&dumpre,&dum);

         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = dumrho;
            B(PRE,i,j) = dumpre;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = l_0*ftheta(t)*(alpha - beta*vr);
         }
      }
      fclose(file);
   }

#endif
}
