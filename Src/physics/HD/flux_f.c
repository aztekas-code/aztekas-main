/**
 * @file /HD/flux_f.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Computes the value of the integration surfaces along
 * all faces of the control volume.
 */

#include"main.h"
    
void Prim2FluxF(double *f, double *v, double *u, gauge_ *local_grid)
{

   double E;
   eos_ eos;
   double rho, p, vx1=0, vx2=0, vx3=0;
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

   EoS(&eos,u,local_grid);

   E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

   f[DEN] = rho * vx1;
   f[ENE] = vx1 * (E + p);
   f[MX1] = rho * vx1 * vx1 + p;
   f[MX2] = rho * vx2 * vx1;
   f[MX3] = rho * vx3 * vx1;

   v[0] = vx1 - eos.cs;
   v[1] = vx1 + eos.cs;
   v[2] = vx1;
}
