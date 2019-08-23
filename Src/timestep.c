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
   int i, j, k, n;
   double dtmin;
   double c, dt, cmax;

   dtmin = 100000.0;

#if DIM == 1

   for(i = gc; i <= Nx1-gc; i++)
   {
      c = sqrt(K*U(1,i) / (U(0,i)));
      dtmin = MIN(dx1/(fabs(U(2,i)) + fabs(c)),dtmin);

      if(U(0,i) == fabs(sqrt(1.0)) || U(0,i) == fabs(1.0/0.0))
      {
         printf("                                          \n");
         printf("NaN value found in calculation.\n");
         exit(EXIT_FAILURE);
      }
   }

#elif DIM == 2

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         c = sqrt(K*U(1,i,j) / (U(0,i,j)));
         dtmin = MIN(dx1/(fabs(U(2,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(3,i,j)) + fabs(c)),dtmin);

         if(U(0,i,j) == fabs(sqrt(1.0)) || U(0,i,j) == fabs(1.0/0.0))
         {
            printf("                                          \n");
            printf("NaN value found in calculation.\n");
            CHECK_NAN = TRUE;
            U = U0;
            PrintValues(&grid.time,&c,&CHECK_NAN);
            exit(EXIT_FAILURE);
         }
      }
   }

#elif DIM == 4

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         c = sqrt(K*U(1,i,j) / (U(0,i,j)));
         dtmin = MIN(dx1/(fabs(U(2,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx1/(fabs(U(4,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(3,i,j)) + fabs(c)),dtmin);
         dtmin = MIN(dx2/(fabs(U(4,i,j)) + fabs(c)),dtmin);

         if(U(0,i,j) == fabs(sqrt(1.0)) || U(0,i,j) == fabs(1.0/0.0))
         {
            printf("                                          \n");
            printf("NaN value found in calculation.\n");
            CHECK_NAN = TRUE;
            U = U0;
            PrintValues(&tprint,&dtprint,&itprint);
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
            dtmin = MIN(dx1/(fabs(U(2,i,j,k)) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U(3,i,j,k)) + fabs(c)),dtmin);
            dtmin = MIN(dx3/(fabs(U(4,i,j,k)) + fabs(c)),dtmin);

            if(U(0,i,j,k) == fabs(sqrt(1.0)) || U(0,i,j,k) == fabs(1.0/0.0))
            {
               printf("                                          \n");
               printf("NaN value found in calculation.\n");
               CHECK_NAN = TRUE;
               U = U0;
               PrintValues(&tprint,&dtprint,&itprint);
               exit(EXIT_FAILURE);
            }
         }
      }
   }
   

#endif

#if PHYSICS == HD
   dt = cou*dtmin;
#elif PHYSICS == RHD
   #if DIM == 1
   dt = cou*MIN(dx1,1000);
   #elif DIM == 2 || DIM == 4
   dt = cou*MIN(dx1,dx2);
   #endif 
#endif

   return dt;
}
