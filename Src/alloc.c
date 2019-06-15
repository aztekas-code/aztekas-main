/**
 * @file alloc.c
 * @author Alejandro Aguayo-Ortiz
 * @brief Essential allocation functions for \a aztekas.
 */

#include"main.h"

void Allocate_Array()
{
#if DIM == 1
   X1  = (double *)malloc((Nx1+1)*sizeof(double));
   X1p = (double *)malloc((Nx1+1)*sizeof(double));
   X1m = (double *)malloc((Nx1+1)*sizeof(double));
   Q   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U   = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(eq+1)*sizeof(double));
#elif DIM == 2 || DIM == 4
   X1  = (double *)malloc((Nx1+1)*sizeof(double));
   X1p = (double *)malloc((Nx1+1)*sizeof(double));
   X1m = (double *)malloc((Nx1+1)*sizeof(double));
   X2  = (double *)malloc((Nx2+1)*sizeof(double));
   X2p = (double *)malloc((Nx2+1)*sizeof(double));
   X2m = (double *)malloc((Nx2+1)*sizeof(double));
   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q   = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q1  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
   Q2  = (double *)malloc((Nx1+1)*(Nx2+1)*(eq+1)*sizeof(double));
#elif DIM == 3
   X1  = (double *)malloc((Nx1+1)*sizeof(double));
   X1p = (double *)malloc((Nx1+1)*sizeof(double));
   X1m = (double *)malloc((Nx1+1)*sizeof(double));
   X2  = (double *)malloc((Nx2+1)*sizeof(double));
   X2p = (double *)malloc((Nx2+1)*sizeof(double));
   X2m = (double *)malloc((Nx2+1)*sizeof(double));
   X3  = (double *)malloc((Nx3+1)*sizeof(double));
   X3p = (double *)malloc((Nx3+1)*sizeof(double));
   X3m = (double *)malloc((Nx3+1)*sizeof(double));
   U   = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
   U1  = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
   U2  = (double *)malloc((Nx1+1)*(Nx2+1)*(Nx3+1)*(eq+1)*sizeof(double));
#endif
}

void New_Size() 
{
   Nx1 = Nx1 + 2*gc;
   Nx2 = Nx2 + 2*gc;
   Nx3 = Nx3 + 2*gc;
}
