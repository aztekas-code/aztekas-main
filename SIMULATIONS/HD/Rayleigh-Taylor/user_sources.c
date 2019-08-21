/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:50:45
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

void User_Source_Terms(double *s, double *u, gauge_ local_grid)
{
   int i;
   double rho, p, vx1=0.0, vx2=0.0, vx3=0.0;

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

   s[0] = 0.0;
   s[1] = -rho*vx2;
   s[2] = 0.0;
   s[3] = -rho;
   s[4] = 0.0;
}
