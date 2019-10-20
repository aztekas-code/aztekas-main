/**
 * @file fluxes.c
 * 
 * @author Alejandro Aguayo-Ortz
 *
 * @brief Numerical fluxes routines.
 *
 * Compute the numerical fluxes at the interfaces for the different
 * directions of the computation
 */

#include"main.h"

void Numerical_Flux_F(double *F, int pm, int *I)
{
   int n, i, j, k;
   double lr, ll;
   double dp[3];
   double dm[3];
   flx_ f;
   gauge_ local_grid;

   local_grid.x[0] = grid.time;

#if DIM == 1

   i = I[0];

   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = 0.0;
   local_grid.x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   local_grid.x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   i = I[0];
   j = I[1];

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

   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = grid.X2[j];
   local_grid.x[3] = grid.X3[k];

#endif

   if(pm == PLUS)
   {
      local_grid.x[1] = grid.X1p[i];

      for(n = 0; n < eq; n++)
      {
         #if DIM == 1
         f.uR[n] = U1m(n,i+1);
         f.uL[n] = U1p(n, i );
         #elif DIM == 2 || DIM == 4
         f.uR[n] = U1m(n,i+1,j);
         f.uL[n] = U1p(n, i ,j);
         #elif DIM == 3
         f.uR[n] = U1m(n,i+1,j,k);
         f.uL[n] = U1p(n, i ,j,k);
         #endif
      }
   }
   if(pm == MINUS)
   {
      local_grid.x[1] = grid.X1m[i];

      for(n = 0; n < eq; n++)
      {
         #if DIM == 1
         f.uR[n] = U1m(n, i );
         f.uL[n] = U1p(n,i-1);
         #elif DIM == 2 || DIM ==4
         f.uR[n] = U1m(n, i ,j);
         f.uL[n] = U1p(n,i-1,j);
         #elif DIM == 3
         f.uR[n] = U1m(n, i ,j,k);
         f.uL[n] = U1p(n,i-1,j,k);
         #endif
      }
   }

#if PHYSICS == RHD
   Get_Metric_Components(&local_grid);
#endif

   Prim2Cons(f.qR,f.uR,&local_grid);
   Prim2Cons(f.qL,f.uL,&local_grid);

   Prim2FluxF(f.fR,dp,f.uR,&local_grid);
   Prim2FluxF(f.fL,dm,f.uL,&local_grid);

   lr = Max(0.0,dp[0],dp[1],dp[2]);
   ll = Max(0.0,dm[0],dm[1],dm[2]);

   f.lR = MAX(lr,ll);

   lr = Min(0.0,dp[0],dp[1],dp[2]);
   ll = Min(0.0,dm[0],dm[1],dm[2]);

   f.lL = MIN(lr,ll);

   Riemann_Solver(F,&f,1);
}

void Numerical_Flux_G(double *F, int pm, int *I)
{
   int n, i, j, k;
   double lr, ll;
   double dp[3];
   double dm[3];
   flx_ f;
   gauge_ local_grid;

   local_grid.x[0] = grid.time;

#if DIM == 2 || DIM == 4

   i = I[0];
   j = I[1];

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

   local_grid.x[1] = grid.X1[i];
   local_grid.x[2] = grid.X2[j];
   local_grid.x[3] = grid.X3[k];

#endif

   if(pm == PLUS)
   {
      local_grid.x[1] = grid.X2p[j];

      for(n = 0; n < eq; n++)
      {
         #if DIM == 2 || DIM == 4
         f.uR[n] = U2m(n,i,j+1);
         f.uL[n] = U2p(n,i, j );
         #elif DIM == 3
         f.uR[n] = U2m(n,i,j+1,k);
         f.uL[n] = U2p(n,i, j ,k);
         #endif
      }
   }
   if(pm == MINUS)
   {
      local_grid.x[1] = grid.X2m[j];

      for(n = 0; n < eq; n++)
      {
         #if DIM == 2 || DIM ==4
         f.uR[n] = U2m(n,i, j );
         f.uL[n] = U2p(n,i,j-1);
         #elif DIM == 3
         f.uR[n] = U2m(n,i, j ,k);
         f.uL[n] = U2p(n,i,j-1,k);
         #endif
      }
   }


#if PHYSICS == RHD
   Get_Metric_Components(&local_grid);
#endif

   Prim2Cons(f.qR,f.uR,&local_grid);
   Prim2Cons(f.qL,f.uL,&local_grid);

   Prim2FluxG(f.fR,dp,f.uR,&local_grid);
   Prim2FluxG(f.fL,dm,f.uL,&local_grid);

   lr = Max(0.0,dp[0],dp[1],dp[2]);
   ll = Max(0.0,dm[0],dm[1],dm[2]);

   f.lR = MAX(lr,ll);

   lr = Min(0.0,dp[0],dp[1],dp[2]);
   ll = Min(0.0,dm[0],dm[1],dm[2]);

   f.lL = MIN(lr,ll);

   Riemann_Solver(F,&f,2);
}
