/*
 * File Name : user_sources.c
 * Description : aztekas user source terms for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 10:01:11
 * Created By :
 */

#include"main.h"

void User_Source_Terms(double *s, double *u, gauge_ local_grid)
{
   double rho, p, vx1=0.0, vx2=0.0, vx3=0.0;
   double r;

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

   r = local_grid.x[1];

   s[RHO] = 0.0;
   s[PRE] = -rho*vx1/(r*r);
   s[VX1] = -rho/(r*r);
   s[VX2] = 0.0;
   s[VX3] = 0.0;
}
