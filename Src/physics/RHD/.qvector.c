#include"main.h"
    
void Prim2Cons(double *q, double *u, gauge_ local_grid)
{
   double n, p, vx1=0, vx2=0, vx3=0;
   double W, h;
   n   = u[0];
   p   = u[1];
   vx1 = u[2];

   W = 1/sqrt((-pow(vx3,2.0))-pow(vx2,2.0)-pow(vx1,2.0)+1);
   h = (K*p+(K-1)*n)/((K-1)*n);

   q[0] = W*n;
   q[1] = (pow(W,2.0)*h - W)*n - p;
   q[2] = pow(W,2.0)*h*n*vx1;
   q[3] = pow(W,2.0)*h*n*vx2;
   q[4] = pow(W,2.0)*h*n*vx3;
}
