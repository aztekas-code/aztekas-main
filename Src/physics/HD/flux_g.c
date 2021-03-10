/**
 * @file /HD/flux_g.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Function that compute Fluxes G along x2 axis 
 * using the primitive variables in HD.
 */

#include"main.h"
    
void Prim2FluxG(double *f, double *v, double *u, gauge_ *local_grid)
{
   double E;
   double rho, p, vx1=0, vx2=0, vx3=0;
   double P[3];
   eos_ eos;
   rho = u[0];
   p   = u[1];

#if DIM == 1
   vx1 = u[2];
#elif DIM == 2
   vx1 = u[2];
   vx2 = u[3];
#elif DIM == 3 || DIM == 4
   vx1 = u[2];
   vx2 = u[3];
   vx3 = u[4];
#endif

   P[0] = rho;
   P[1] = p;
   P[2] = 0.0;
   EoS(&eos,P,local_grid);

   E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

   f[DEN] = rho * vx2;
   f[ENE] = vx2 * (E + p);
   f[MX1] = rho * vx1 * vx2;
   f[MX2] = rho * vx2 * vx2 + p;
   f[MX3] = rho * vx3 * vx2;

   v[0] = vx2 - eos.cs;
   v[1] = vx2 + eos.cs;
   v[2] = vx2;
}
