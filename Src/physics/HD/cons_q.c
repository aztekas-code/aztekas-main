/**
 * @file /HD/cons_q.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Function that converts Primitives to Conservative variables in HD.
 */

#include"main.h"
    
void Prim2Cons(double *q, double *u, gauge_ *local_grid)
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

   q[DEN] = rho;
   q[ENE] = E;
   q[MX1] = rho*vx1;
   q[MX2] = rho*vx2;
   q[MX3] = rho*vx3;
}
