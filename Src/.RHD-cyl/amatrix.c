#include<stdio.h>
#include<math.h>
#include"../Headers/vector.h"
#include"../Headers/main.h"

int funct_A(double *a, double *uu)
{
   int i, j;
   double n, p, u=0, v=0, w=0;
   double R, W, h;
   double dWu, dWv, dWw;
   double dhn, dhp;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}

   R = sqrt(pow(x2,2.0)+pow(x1,2.0));
   W = x1/sqrt(-pow(w,2.0)-pow(x1,2.0)*pow(v,2.0)-pow(x1,2.0)*pow(u,2.0)+pow(x1,2.0));
   h = (K*p+(K-1)*n)/((K-1)*n);

   dWu = u*pow(W,3.0);
   dWv = v*pow(W,3.0);
   dWw = w*pow(W,3.0)/pow(x1,2.0);
   dhn = -K*p/((K-1)*pow(n,2.0));
   dhp = K/((K-1)*n);

   a[0*eq + 0] = ((dhp*n*w*W-2*h*w)*dWw+(dhp*n*v*W-2*h*v)*dWv+(dhp*n*u*W-2*h*u)*dWu+h*dhp*n*pow(W,3.0)-h*W)/((dhn*n-h)*w*W*dWw+(dhn*n-h)*v*W*dWv+(dhn*n-h)*u*W*dWu+h*dhp*n*pow(W,4.0)-h*pow(W,2.0));
   a[0*eq + 1] = (dhp*n*w*dWw+dhp*n*v*dWv+dhp*n*u*dWu)/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[0*eq + 2] = -(dhp*n*pow(W,2.0)-1)*dWu/((dhn*n-h)*w*pow(W,2.0)*dWw+(dhn*n-h)*v*pow(W,2.0)*dWv+(dhn*n-h)*u*pow(W,2.0)*dWu+h*dhp*n*pow(W,5.0)-h*pow(W,3.0));
   a[0*eq + 3] = -(dhp*n*pow(W,2.0)-1)*dWv/((dhn*n-h)*w*pow(W,2.0)*dWw+(dhn*n-h)*v*pow(W,2.0)*dWv+(dhn*n-h)*u*pow(W,2.0)*dWu+h*dhp*n*pow(W,5.0)-h*pow(W,3.0));
   a[0*eq + 4] = -(dhp*n*pow(W,2.0)-1)*dWw/((dhn*n-h)*w*pow(W,2.0)*dWw+(dhn*n-h)*v*pow(W,2.0)*dWv+(dhn*n-h)*u*pow(W,2.0)*dWu+h*dhp*n*pow(W,5.0)-h*pow(W,3.0));
   a[1*eq + 0] = -((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+(h*dhn*n+pow(h,2.0))*pow(W,2.0)-h*W)/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[1*eq + 1] = -((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu-h*W)/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[1*eq + 2] = (dhn*n-h)*dWu/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[1*eq + 3] = (dhn*n-h)*dWv/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[1*eq + 4] = (dhn*n-h)*dWw/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[2*eq + 0] = -(dhp*n*u*W+(-dhn*n-h)*u)/((dhn*pow(n,2.0)-h*n)*w*dWw+(dhn*pow(n,2.0)-h*n)*v*dWv+(dhn*pow(n,2.0)-h*n)*u*dWu+h*dhp*pow(n,2.0)*pow(W,3.0)-h*n*W);
   a[2*eq + 1] = -dhp*u*W/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[2*eq + 2] = ((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+h*dhp*n*pow(W,3.0)-h*W)/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[2*eq + 3] = -(dhn*n-h)*u*dWv/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[2*eq + 4] = -(dhn*n-h)*u*dWw/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[3*eq + 0] = -(dhp*n*v*W+(-dhn*n-h)*v)/((dhn*pow(n,2.0)-h*n)*w*dWw+(dhn*pow(n,2.0)-h*n)*v*dWv+(dhn*pow(n,2.0)-h*n)*u*dWu+h*dhp*pow(n,2.0)*pow(W,3.0)-h*n*W);
   a[3*eq + 1] = -dhp*v*W/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[3*eq + 2] = -(dhn*n-h)*v*dWu/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[3*eq + 3] = ((dhn*n-h)*w*dWw+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W)/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[3*eq + 4] = -(dhn*n-h)*v*dWw/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[4*eq + 0] = -(dhp*n*w*W+(-dhn*n-h)*w)/((dhn*pow(n,2.0)-h*n)*w*dWw+(dhn*pow(n,2.0)-h*n)*v*dWv+(dhn*pow(n,2.0)-h*n)*u*dWu+h*dhp*pow(n,2.0)*pow(W,3.0)-h*n*W);
   a[4*eq + 1] = -dhp*w*W/((dhn*n-h)*w*dWw+(dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W);
   a[4*eq + 2] = -(dhn*n-h)*w*dWu/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[4*eq + 3] = -(dhn*n-h)*w*dWv/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));
   a[4*eq + 4] = ((dhn*n-h)*v*dWv+(dhn*n-h)*u*dWu+h*dhp*n*pow(W,3.0)-h*W)/((h*dhn*pow(n,2.0)-pow(h,2.0)*n)*w*pow(W,2.0)*dWw+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*v*pow(W,2.0)*dWv+(h*dhn*pow(n,2.0)-pow(h,2.0)*n)*u*pow(W,2.0)*dWu+pow(h,2.0)*dhp*pow(n,2.0)*pow(W,5.0)-pow(h,2.0)*n*pow(W,3.0));

   return 0;
}
