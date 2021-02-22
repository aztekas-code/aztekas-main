/**
 * @file /HD/flux_h.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Function that compute Fluxes H along x3 axis 
 * using the primitive variables in HD.
 */

#include"main.h"
    
//void Prim2FluxH(double *f, double *v, double *u, double *x)
void Prim2FluxH(double *f, double *v, double *u, gauge_ *local_grid)
{
   double E, cs;
   eos_ eos;
   double rho, p, vx1=0, vx2=0, vx3=0;
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

   EoS(&eos,u,local_grid);

   E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

   f[0] = rho * vx3;
   f[1] = vx3 * (E + p);
   f[2] = rho * vx1 * vx3;
   f[3] = rho * vx2 * vx3;
   f[4] = rho * vx3 * vx3 + p;

   v[0] = vx3 - cs;
   v[1] = vx3 + cs;
   v[2] = vx3;
}
