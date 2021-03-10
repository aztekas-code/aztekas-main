/**
 * @file /HD/flux_h.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Function that compute Fluxes H along x3 axis 
 * using the primitive variables in HD.
 */

#include"main.h"
    
void Prim2FluxH(double *f, double *v, double *u, gauge_ *local_grid)
{
   double E;
   double rho, p, vx1=0, vx2=0, vx3=0;
   double P[3];
   eos_ eos;
   rho = u[RHO];
   p   = u[PRE];

#if DIM == 1
   vx1 = u[VX1];
#elif DIM == 2
   vx1 = u[VX1];
   vx2 = u[VX2];
#elif DIM == 3 || DIM == 4
   vx1 = u[VX1];
   vx2 = u[VX2];
   vx3 = u[VX3];
#endif

   P[0] = rho;
   P[1] = p;
   P[2] = 0.0;
   EoS(&eos,P,local_grid);

   E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

   f[DEN] = rho * vx3;
   f[ENE] = vx3 * (E + p);
   f[MX1] = rho * vx1 * vx3;
   f[MX2] = rho * vx2 * vx3;
   f[MX3] = rho * vx3 * vx3 + p;

   v[0] = vx3 - eos.cs;
   v[1] = vx3 + eos.cs;
   v[2] = vx3;
}
