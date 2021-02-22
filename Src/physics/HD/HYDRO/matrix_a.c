/**
 * @file /HD/matrix_a.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Computes the matrix \f$ A := (\partial {Q}/\partial
 * {U}) \f$, needed for the PVRS method.
 */

#include"main.h"
    
void Matrix_A(double *a, double *u, gauge_ local_grid)
{
   double rho, p, vx1, vx2, vx3;

   // Density and Pressure
   rho = u[RHO];
   p   = u[PRE];

   // Covariant components of the 3-velocity
#if DIM == 1
   vx1 = u[VX1];
   vx2 = 0.0;
   vx3 = 0.0;
#elif DIM == 2
   vx1 = u[VX1];
   vx2 = u[VX2];
   vx3 = 0.0;
#elif DIM == 3 || DIM == 4
   vx1 = u[VX1];
   vx2 = u[VX2];
   vx3 = u[VX3];
#endif

#if DIM == 1

   a(0,0) = 1.0;
   a(0,1) = 0.0;
   a(0,2) = 0.0;

   a(1,0) = ((K-1)*pow(vx3,2.0)+(K-1)*pow(vx2,2.0)+(K-1)*pow(vx1,2.0))/2;
   a(1,1) = K-1;
   a(1,2) = (1-K)*vx1;

   a(2,0) = -vx1/rho;
   a(2,1) = 0;
   a(2,2) = 1.0/rho;

#elif DIM == 2 

   a(0,0) = 1.0;
   a(0,1) = 0.0;
   a(0,2) = 0.0;
   a(0,3) = 0.0;

   a(1,0) = ((K-1)*pow(vx3,2.0)+(K-1)*pow(vx2,2.0)+(K-1)*pow(vx1,2.0))/2;
   a(1,0) = K-1;
   a(1,2) = (1-K)*vx1;
   a(1,3) = (1-K)*vx2;

   a(2,0) = -vx1/rho;
   a(2,1) = 0;
   a(2,2) = 1.0/rho;
   a(2,3) = 0;

   a(3,0) = -vx2/rho;
   a(3,1) = 0;
   a(3,2) = 0;
   a(3,3) = 1.0/rho;

#elif DIM == 3 || DIM == 4

   a(0,0) = 1.0;
   a(0,1) = 0.0;
   a(0,2) = 0.0;
   a(0,3) = 0.0;
   a(0,4) = 0.0;

   a(1,0) = ((K-1)*pow(vx3,2.0)+(K-1)*pow(vx2,2.0)+(K-1)*pow(vx1,2.0))/2;
   a(1,0) = K-1;
   a(1,2) = (1-K)*vx1;
   a(1,3) = (1-K)*vx2;
   a(1,4) = (1-K)*vx3;

   a(2,0) = -vx1/rho;
   a(2,1) = 0;
   a(2,2) = 1.0/rho;
   a(2,3) = 0;
   a(2,4) = 0;

   a(3,0) = -vx2/rho;
   a(3,1) = 0;
   a(3,2) = 0;
   a(3,3) = 1.0/rho;
   a(3,4) = 0;

   a(4,0) = -vx3/rho;
   a(4,1) = 0;
   a(4,2) = 0;
   a(4,3) = 0;
   a(4,4) = 1.0/rho;

#endif
}
