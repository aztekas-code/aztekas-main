/**
 * @file /grid/mesh.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Cartesian-like mesh grid construction.
 */

#include"main.h"

int Mesh()
{
   int i, j, k;
   double M, a, dr, Delta;

#if DIM == 1
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
#elif DIM == 2 || DIM == 4
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
#elif DIM == 3
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
   dx3 = (x3max - x3min)/((double)Nx3-2*gc);
#endif

#if DIM == 1
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i] = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);

      #if GRID == LOGMESH
      grid.X1[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2*gc)) - 1;
      grid.X1p[i] = x1min + exp(log((x1max - x1min + 1.0))*(i+0.5-gc)/(Nx1-2*gc)) - 1;
      grid.X1m[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-0.5-gc)/(Nx1-2*gc)) - 1;
      #endif
   }

#if PHYSICS == HD || PHYSICS == RHD
   Surface_Volume();
#endif

#elif DIM == 2  || DIM == 4
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i]  = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);

      #if GRID == LOGMESH
      grid.X1[i]  = x1min + lfac*exp(log((x1max - x1min + lfac)/lfac)*(i-gc)/(Nx1-2*gc))     - lfac;
      grid.X1p[i] = x1min + lfac*exp(log((x1max - x1min + lfac)/lfac)*(i+0.5-gc)/(Nx1-2*gc)) - lfac;
      grid.X1m[i] = x1min + lfac*exp(log((x1max - x1min + lfac)/lfac)*(i-0.5-gc)/(Nx1-2*gc)) - lfac;
//      grid.X1[i]  = x1min + lfac*exp(log((5*rs(K,cs) - x1min + lfac)/lfac)*(i-gc)/(100))     - lfac;
//      grid.X1p[i] = x1min + lfac*exp(log((5*rs(K,cs) - x1min + lfac)/lfac)*(i+0.5-gc)/(100)) - lfac;
//      grid.X1m[i] = x1min + lfac*exp(log((5*rs(K,cs) - x1min + lfac)/lfac)*(i-0.5-gc)/(100)) - lfac;
//      grid.X1[i]  = x1min + lfac*exp(log((1*rB - x1min + lfac)/lfac)*(i-gc)/(100))     - lfac;
//      grid.X1p[i] = x1min + lfac*exp(log((1*rB - x1min + lfac)/lfac)*(i+0.5-gc)/(100)) - lfac;
//      grid.X1m[i] = x1min + lfac*exp(log((1*rB - x1min + lfac)/lfac)*(i-0.5-gc)/(100)) - lfac;
      #endif
   }

   for(j = 0; j <= Nx2; j++)
   {
      grid.X2[j]  = x2min + (j-gc)*(dx2);
      grid.X2p[j] = x2min + (j+0.5-gc)*(dx2);
      grid.X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }
   
#if PHYSICS == HD || PHYSICS == RHD
   Surface_Volume();
#endif

#elif DIM == 3 
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i]  = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);
   }

   for(j = 0; j <= Nx2; j++)
   {
      grid.X2[j]  = x2min + (j-gc)*(dx2);
      grid.X2p[j] = x2min + (j+0.5-gc)*(dx2);
      grid.X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }

   for(k = 0; k <= Nx3; k++)
   {
      grid.X3[k]  = x3min + (k-gc)*(dx3);
      grid.X3p[k] = x3min + (k+0.5-gc)*(dx3);
      grid.X3m[k] = x3min + (k-0.5-gc)*(dx3);
   }

#endif
   return 0;
}
