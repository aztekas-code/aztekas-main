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
   int i, j, k;
   double dtmin;
   double c, dt, cmax;
   double r;

   dtmin = 100000;

   if(DIM == 1)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         c = sqrt(K*U[c1(1,i)] / (U[c1(0,i)]));
         dtmin = MIN(dx1/(fabs(U[c1(2,i)]) + fabs(c)),dtmin);
      }
   }
   else if(DIM == 2)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            c = sqrt(K*U[c2(1,i,j)] / (U[c2(0,i,j)]));
            dtmin = MIN(dx1/(fabs(U[c2(2,i,j)]) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U[c2(3,i,j)]) + fabs(c)),dtmin);
         }
      }
   }
   else if(DIM == 4)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            c = sqrt(K*U[c2(1,i,j)] / (U[c2(0,i,j)]));
            dtmin = MIN(dx1/(fabs(U[c2(2,i,j)]) + fabs(c)),dtmin);
            dtmin = MIN(dx1/(fabs(U[c2(4,i,j)]) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U[c2(3,i,j)]) + fabs(c)),dtmin);
            dtmin = MIN(dx2/(fabs(U[c2(4,i,j)]) + fabs(c)),dtmin);
         }
      }
   }
   else if(DIM == 3)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            for(k = gc; k <= Nx3-gc; k++)
            {
               c = sqrt(K*U[c3(1,i,j,k)] / (U[c3(0,i,j,k)]));
               dtmin = MIN(dx1/(fabs(U[c3(2,i,j,k)]) + fabs(c)),dtmin);
               dtmin = MIN(dx2/(fabs(U[c3(3,i,j,k)]) + fabs(c)),dtmin);
               dtmin = MIN(dx3/(fabs(U[c3(4,i,j,k)]) + fabs(c)),dtmin);
            }
         }
      }
   }

#if PHYSICS == 1 //HD
   dt = cou*dtmin;
#elif PHYSICS == 2 //RHD
   #if DIM == 1
   dt = cou*MIN(dx1,1000);
   #elif DIM == 2 || DIM == 4
   dt = cou*MIN(dx1,dx2);
   #endif 
#endif

   return dt;
}
