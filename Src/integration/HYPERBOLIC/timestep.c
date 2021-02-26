/**
 * @file timestep.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Time-step calculation.
 */

//Do not erase any of these libraries//
#include"main.h"

double TimeStep()
{
   double dtmin, Dx1, Dx2;
   double c, dt, cmax;

   dtmin = 100000.0;

#if DIM == 1

#ifdef _OPENMP
   #pragma omp parallel shared(U,grid)
   #pragma omp for private(c,Dx1) reduction(min : dtmin)
#endif
   for(int i = gc; i <= Nx1-gc; i++)
   {
      Dx1 = grid.X1p[i] - grid.X1m[i];

      #if PHYSICS == HD
      c = sqrt(K*U(PRE,i) / (U(RHO,i)));

      dtmin = MIN(Dx1/(fabs(U(VX1,i)) + fabs(c)),dtmin);
      #elif PHYSICS == RHD
      dtmin = MIN(Dx1,dtmin);
      #endif

      if(U(0,i) != U(0,i) || U(1,i) != U(1,i))
      {
         printf("                                          \n");
         printf("NaN value found in calculation at (%d).\n",i);
         CHECK_NAN = TRUE;
         U = U0;
         Print_Time_Values(&grid.time,&c,&CHECK_NAN);
         exit(EXIT_FAILURE);
      }
   }


#elif DIM == 2

#ifdef _OPENMP
   #pragma omp parallel shared(U,grid) 
   #pragma omp for private(c,Dx1,Dx2) reduction(min : dtmin) collapse(2)
#endif
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         Dx1 = grid.X1p[i] - grid.X1m[i];
         Dx2 = grid.X2p[j] - grid.X2m[j];
         #if COORDINATES == SPHERICAL
         if(grid.X1[i] != 0.0)
         {
            Dx2 = grid.X1[i]*Dx2;
         }
         else
         {
            Dx2 = Dx2;
         }
         #endif

         #if PHYSICS == HD
         c = sqrt(K*U(PRE,i,j) / (U(RHO,i,j)));

         dtmin = MIN(Dx1/(fabs(U(VX1,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(Dx2/(fabs(U(VX2,i,j)) + fabs(c)),dtmin);
         #elif PHYSICS == RHD
         dtmin = MIN(Dx1,dtmin);
         dtmin = MIN(Dx2,dtmin);
         #endif
         
         if(U(0,i,j) != U(0,i,j) || U(1,i,j) != U(1,i,j))
         {
            printf("                                          \n");
            printf("NaN value found in calculation at (%d,%d) = (%.2f,%.2f).\n",i,j,grid.X1[i],grid.X2[j]);
            CHECK_NAN = TRUE;
            U = U0;
            Print_Time_Values(&grid.time,&c,&CHECK_NAN);
            exit(EXIT_FAILURE);
         }
      }
   }

#elif DIM == 4

#ifdef _OPENMP
   #pragma omp parallel shared(U,grid)
   #pragma omp for private(c,Dx1,Dx2) reduction(min : dtmin) collapse(2)
#endif
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         Dx1 = grid.X1p[i] - grid.X1m[i];
         Dx2 = grid.X2p[j] - grid.X2m[j];
         #if COORDINATES == SPHERICAL
         Dx2 = grid.X1[i]*Dx2;
         #endif

         #if PHYSICS == HD
         c = sqrt(K*U(PRE,i,j) / (U(RHO,i,j)));

         dtmin = MIN(Dx1/(fabs(U(VX1,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(Dx1/(fabs(U(VX3,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(Dx2/(fabs(U(VX2,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(Dx2/(fabs(U(VX3,i,j)) + fabs(c)),dtmin);
         #elif PHYSICS == RHD
         dtmin = MIN(Dx1,dtmin);
         dtmin = MIN(Dx2,dtmin);
         #endif

         if(U(0,i,j) != U(0,i,j) || U(1,i,j) != U(1,i,j))
         {
            printf("                                          \n");
            printf("NaN value found in calculation at (%d,%d).\n",i,j);
            CHECK_NAN = TRUE;
            U = U0;
            Print_Time_Values(&grid.time,&c,&CHECK_NAN);
            exit(EXIT_FAILURE);
         }
      }
   }

#elif DIM == 3

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            c = sqrt(K*U(1,i,j,k) / (U(0,i,j,k)));

            Dx1 = grid.X1p[i] - grid.X1m[i];
            Dx2 = grid.X2p[j] - grid.X2m[j];
            Dx3 = grid.X3p[k] - grid.X3m[k];
            dtmin = MIN(dx1/(fabs(U(2,i,j,k)) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U(3,i,j,k)) + fabs(c)),dtmin);
            dtmin = MIN(dx3/(fabs(U(4,i,j,k)) + fabs(c)),dtmin);

//            #pragma omp critical
            if(U(RHO,i,j,k) == fabs(1.0/0.0))
            {
               printf("                                          \n");
               printf("NaN value found in calculation at (%d,%d,%d).\n",i,j,k);
               CHECK_NAN = TRUE;
               U = U0;
               Print_Time_Values(&tprint,&dtprint,&itprint);
               exit(EXIT_FAILURE);
            }
         }
      }
   }
   

#endif

   dt = cou*dtmin;

   return dt;
}
