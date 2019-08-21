/**
 * @file alloc.c
 * @author Alejandro Aguayo-Ortiz
 * @brief Essential allocation functions for \a aztekas.
 */

#include"main.h"

void Allocate_Array()
{
#if DIM == 1
   grid.X1 = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));

   grid.S1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.S1m = (double *)malloc((Nx1+1)*sizeof(double));

   Q   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
#elif DIM == 2 || DIM == 4
   grid.X1  = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X2  = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2p = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2m = (double *)malloc((Nx2+1)*sizeof(double));

   grid.S1p = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S1m = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S2p = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));
   grid.S2m = (double *)malloc((Nx1+1)*(Nx2+1)*sizeof(double));

   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
#elif DIM == 3
   grid.X1  = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1p = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X1m = (double *)malloc((Nx1+1)*sizeof(double));
   grid.X2  = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2p = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X2m = (double *)malloc((Nx2+1)*sizeof(double));
   grid.X3  = (double *)malloc((Nx3+1)*sizeof(double));
   grid.X3p = (double *)malloc((Nx3+1)*sizeof(double));
   grid.X3m = (double *)malloc((Nx3+1)*sizeof(double));
   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
#endif
}

/*
void Allocate_Metric()
{
   int i, j, k;
#if DIM == 1
   for(i = 0; i < DIM + 1; i++)
   {
      for(j = 0; j < DIM + 1; j++)
      {
         gmetric.g = (double *)malloc((Nx1+1)*sizeof(double));
      }
   }

   for(i = 0; i < DIM; i++)
   {
      for(j = 0; j < DIM; j++)
      {

      }
   }


}
*/

void New_Size()
{
   Nx1 = Nx1 + 2*gc;
   Nx2 = Nx2 + 2*gc;
   Nx3 = Nx3 + 2*gc;
}
