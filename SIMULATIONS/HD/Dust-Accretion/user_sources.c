/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:13:33
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

void User_Source_Terms(double *a, double *u, double *x)
{
   double n, p, vx1=0.0, vx2=0.0, vx3=0.0;
   n = u[0];
   p = u[1];
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

   a[0] = 0.0;
   a[1] = -n*vx1/(x[1]*x[1]);
   a[2] = -n/(x[1]*x[1]);
   a[3] = 0.0;
   a[4] = 0.0;
}
