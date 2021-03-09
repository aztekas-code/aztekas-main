/**
 * @file /fluxes/sources.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Source terms calculation.
 *
 * Computes the default and user defined source terms.
 */

//Do not erase any of these libraries//
#include"main.h"

void Prim2Sources(double *s, int *I)
{
   int i, j, k;
   double u[10];
   double default_S[10];
#if GRAVITY != NONE
   double grav_S[10];
#endif
#if USER_SOURCE_TERMS == TRUE
   double user_S[10];
#endif
   gauge_ local_grid;

   local_grid.x[0] = grid.time;

#if DIM == 1

   i = I[0];

   local_grid.I[0] = i;
   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = 0.0;
   local_grid.x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   local_grid.x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   i = I[0];
   j = I[1];

   local_grid.I[0] = i;
   local_grid.I[1] = j;
   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = grid.X2[j];
   local_grid.x[3] = 0.0;
   #if POLAR == TRUE
   local_grid.x[2] = M_PI_2;
   #endif

#elif DIM == 3

   i = I[0];
   j = I[1];
   k = I[2];

   local_grid.I[0] = i;
   local_grid.I[1] = j;
   local_grid.I[2] = k;
   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = grid.X2[j];
   local_grid.x[3] = grid.X3[k];

#endif

#if PHYSICS == RHD
   Get_Metric_Components(&local_grid);
#endif

   for(int n = 0; n < eq; n++)
   {
   #if DIM == 1
      u[n] = U(n,i);
   #elif DIM == 2 || DIM == 4
      u[n] = U(n,i,j);
   #elif DIM == 3
      u[n] = U(n,i,j,k);
   #endif
   }

   // Geometric source terms
   Source_Terms(default_S,u,&local_grid);

#if USER_SOURCE_TERMS == TRUE
   User_Source_Terms(user_S,u,&local_grid);
#endif

   for(int n = 0; n < eq; n++)
   {
      s[n] = default_S[n];

   #if USER_SOURCE_TERMS == TRUE
      s[n] += user_S[n];
   #endif
   }

#if INTEGRATION == PVRS
   double matrix_A[(eq+1)*(eq+1)];
   Matrix_A(v->A,u,&local_grid);
#endif

}
