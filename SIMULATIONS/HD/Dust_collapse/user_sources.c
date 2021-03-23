/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:13:33
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

void User_Source_Terms(double *s, double *u, gauge_ *local_grid)
{
   double rho, p, vx1=0.0, vx2=0.0, vx3=0.0;
   double r, mass_sg;
   int iR;

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

   r = local_grid->x[1];

   s[RHO] = 0.0;
   s[VX2] = 0.0;
   s[VX3] = 0.0;

#if selfgravity
   iR = local_grid->I[0];
   mass_sg = mass[iR];
   s[PRE] = -rho*(mass_sg)*vx1/(r*r);
   s[VX1] = -rho*(mass_sg)/(r*r);
#else
   s[PRE] = -rho*vx1/(r*r);
   s[VX1] = -rho/(r*r);   
#endif


}
