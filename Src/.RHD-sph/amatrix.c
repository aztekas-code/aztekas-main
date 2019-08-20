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

   R = x1;
   W = (x1*(sin(x2)))/sqrt(-((x1+2*MM)*pow(w,2.0)+(x1+2*MM)*pow(sin(x2),2.0)*pow(v,2.0)+pow(x1,3.0)*pow(sin(x2),2.0)*pow(u,2.0)+((-pow(x1,3.0))-2*MM*pow(x1,2.0))*pow(sin(x2),2.0))/(x1+2*MM));
   h = (K*p+(K-1)*n)/((K-1)*n);

   dWu = (pow(W,3.0)*x1*u)/(x1+2*MM);
   dWv = (pow(W,3.0)*v)/pow(x1,2.0);
   dWw = (pow(W,3.0)*w)/(pow(x1,2.0)*pow(sin(x2),2.0));
   dhn = -(K*p)/((K-1)*pow(n,2.0));
   dhp = K/((K-1)*n);

   a[0*eq + 0] = ((W*dWw*dhp*n-2*dWw*h)*w+(W*dWv*dhp*n-2*dWv*h)*v+(W*dWu*dhp*n-2*dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h)/((W*dWw*dhn*n-W*dWw*h)*w+(W*dWv*dhn*n-W*dWv*h)*v+(W*dWu*dhn*n-W*dWu*h)*u+pow(W,4.0)*h*dhp*n-pow(W,2.0)*h);
   a[0*eq + 1] = (dWw*dhp*n*w+dWv*dhp*n*v+dWu*dhp*n*u)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[0*eq + 2] = -(pow(W,2.0)*dWu*dhp*n-dWu)/((pow(W,2.0)*dWw*dhn*n-pow(W,2.0)*dWw*h)*w+(pow(W,2.0)*dWv*dhn*n-pow(W,2.0)*dWv*h)*v+(pow(W,2.0)*dWu*dhn*n-pow(W,2.0)*dWu*h)*u+pow(W,5.0)*h*dhp*n-pow(W,3.0)*h);
   a[0*eq + 3] = -(pow(W,2.0)*dWv*dhp*n-dWv)/((pow(W,2.0)*dWw*dhn*n-pow(W,2.0)*dWw*h)*w+(pow(W,2.0)*dWv*dhn*n-pow(W,2.0)*dWv*h)*v+(pow(W,2.0)*dWu*dhn*n-pow(W,2.0)*dWu*h)*u+pow(W,5.0)*h*dhp*n-pow(W,3.0)*h);
   a[0*eq + 4] = -(pow(W,2.0)*dWw*dhp*n-dWw)/((pow(W,2.0)*dWw*dhn*n-pow(W,2.0)*dWw*h)*w+(pow(W,2.0)*dWv*dhn*n-pow(W,2.0)*dWv*h)*v+(pow(W,2.0)*dWu*dhn*n-pow(W,2.0)*dWu*h)*u+pow(W,5.0)*h*dhp*n-pow(W,3.0)*h);
   a[1*eq + 0] = -((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,2.0)*h*dhn*n+pow(W,2.0)*pow(h,2.0)-W*h)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[1*eq + 1] = -((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u-W*h)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[1*eq + 2] = (dWu*dhn*n-dWu*h)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[1*eq + 3] = (dWv*dhn*n-dWv*h)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[1*eq + 4] = (dWw*dhn*n-dWw*h)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[2*eq + 0] = -(((W*dhp-dhn)*n-h)*u)/((dWw*dhn*pow(n,2.0)-dWw*h*n)*w+(dWv*dhn*pow(n,2.0)-dWv*h*n)*v+(dWu*dhn*pow(n,2.0)-dWu*h*n)*u+pow(W,3.0)*h*dhp*pow(n,2.0)-W*h*n);
   a[2*eq + 1] = -(W*dhp*u)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[2*eq + 2] = ((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+pow(W,3.0)*h*dhp*n-W*h)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[2*eq + 3] = -((dWv*dhn*n-dWv*h)*u)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[2*eq + 4] = -((dWw*dhn*n-dWw*h)*u)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[3*eq + 0] = -(((W*dhp-dhn)*n-h)*v)/((dWw*dhn*pow(n,2.0)-dWw*h*n)*w+(dWv*dhn*pow(n,2.0)-dWv*h*n)*v+(dWu*dhn*pow(n,2.0)-dWu*h*n)*u+pow(W,3.0)*h*dhp*pow(n,2.0)-W*h*n);
   a[3*eq + 1] = -(W*dhp*v)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[3*eq + 2] = -((dWu*dhn*n-dWu*h)*v)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[3*eq + 3] = ((dWw*dhn*n-dWw*h)*w+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[3*eq + 4] = -((dWw*dhn*n-dWw*h)*v)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[4*eq + 0] = -(((W*dhp-dhn)*n-h)*w)/((dWw*dhn*pow(n,2.0)-dWw*h*n)*w+(dWv*dhn*pow(n,2.0)-dWv*h*n)*v+(dWu*dhn*pow(n,2.0)-dWu*h*n)*u+pow(W,3.0)*h*dhp*pow(n,2.0)-W*h*n);
   a[4*eq + 1] = -(W*dhp*w)/((dWw*dhn*n-dWw*h)*w+(dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h);
   a[4*eq + 2] = -((dWu*dhn*n-dWu*h)*w)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[4*eq + 3] = -((dWv*dhn*n-dWv*h)*w)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);
   a[4*eq + 4] = ((dWv*dhn*n-dWv*h)*v+(dWu*dhn*n-dWu*h)*u+pow(W,3.0)*h*dhp*n-W*h)/((pow(W,2.0)*dWw*h*dhn*pow(n,2.0)-pow(W,2.0)*dWw*pow(h,2.0)*n)*w+(pow(W,2.0)*dWv*h*dhn*pow(n,2.0)-pow(W,2.0)*dWv*pow(h,2.0)*n)*v+(pow(W,2.0)*dWu*h*dhn*pow(n,2.0)-pow(W,2.0)*dWu*pow(h,2.0)*n)*u+pow(W,5.0)*pow(h,2.0)*dhp*pow(n,2.0)-pow(W,3.0)*pow(h,2.0)*n);

   return 0;
}
