/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"./Headers/main.h"
#include"./Headers/matrix.h"
#include"./Headers/vector.h"

int RK1D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i;
   int I[3];
   double Dx1 = dx1;
   double Dt  = dt;
   double L[eq+1], F[eq+1];
   double UU[eq+1];
   double g[6], g1p[6], g1m[6];
   double g2p[6], g2m[6];
   vec_ v;
   lim_ l;

   #pragma omp parallel shared(u,u1,Dt,Dx1,order,r) private(n,L,F,v,l,i) num_threads(1)
   {
      switch(order)
      {
         case 1:
            #pragma omp for
            for(i = 4; i <= Nx1-4; i++)
            {
               //#pragma omp atomic read
               I[0] = i;
               r = 'C';

               RECONST1D(u,r,&l,I);
               FLUX1D(&v,&l,I);

               gauge(g,X1[i],0,0);
               gauge(g1p,X1p[i],0,0);
               gauge(g1m,X1m[i],0,0);

               for(n = 0; n < eq; n++)
               {
                  F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (v.Fp[n] - v.Fm[n])/(g[1]*Dx1) - \
                  g[0]*v.S[n];
                  //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  g[0]*v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  v.S[n] +\
                  v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) - \
                  v.S[n];
                  //F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  (2.0 - MM*pow(g[0],2.0)/X1[i])*(1.0/X1[i])*v.F[n];
               }

               for(n = 0; n < eq; n++)
               {
                  q1[c1(n,i)] = q[c1(n,i)] - (Dt)*(F[n]);
               }
            }
         break;

         case 2:
            #pragma omp for
            for(i = 4; i <= Nx1-4; i++)
            {
               //#pragma omp atomic read
               I[0] = i;
               r = 'C';

               RECONST1D(u1,r,&l,I);
               FLUX1D(&v,&l,I);

               gauge(g,X1[i],0,0);
               gauge(g1p,X1p[i],0,0);
               gauge(g1m,X1m[i],0,0);

               for(n = 0; n < eq; n++)
               {
                  F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  g[0]*v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
                  v.S[n] +\
                  v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1);
                  //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) - \
                  v.S[n];
                  //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) - \
                  g[0]*v.S[n] +\
                  (2.0 - MM*pow(g[0],2.0)/X1[i])*(g[0]/X1[i])*v.F[n];
               }

               for(n = 0; n < eq; n++)
               {
                  q2[c1(n,i)] = 0.5*(q1[c1(n,i)] + q[c1(n,i)] - (Dt)*F[n]);
               }
            }

         break;
      }
   }

   return 0;
}

int RK2D(double *u, double *u1, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i, j;
   int I[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dt  = dt;
   double L[eq+1], F[eq+1];
   double g[6], g1p[6], g1m[6];
   double g2p[6], g2m[6];
   vec_ v;
   lim_ l;

   #pragma omp parallel shared(u,u1,Dt,Dx1,order,r) private(n,L,F,v,i,j) num_threads(1)
   {
      switch(order)
      {
         case 1:
            #pragma omp for
            for(i = 4; i <= Nx1-4; i++)
            {
               for(j = 4; j <= Nx2-4; j++)
               {
                  I[0] = i;
                  I[1] = j;
                  r = 'G';

                  RECONST2D(u,r,&l,I);
                  FLUX2D(&v,&l,I);

                  gauge(g,X1[i],X2[j],0);
                  gauge(g1p,X1p[i],X2[j],0);
                  gauge(g1m,X1m[i],X2[j],0);
                  gauge(g2p,X1[i],X2p[j],0);
                  gauge(g2m,X1[i],X2m[j],0);

                  for(n = 0; n < eq; n++)
                  {
                     //F[n] = g[1]*(g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     g[1]*(g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[2]*v.S[n] +\
                     (g1p[1] - g1m[1])/(Dx1)*v.F[n] +\
                     (g2p[1] - g2m[1])/(Dx2)*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(g[0]*Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(g[0]*Dx2) -\
                     v.S[n] +\
                     (g1p[1] - g1m[1])/(g[1]*Dx1)*v.F[n] +\
                     (g2p[1] - g2m[1])/(g[1]*Dx2)*v.G[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[2]*Dx2) - v.S[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(Dx2) - g[2]*v.S[n];
                     //F[n] = g[2]*(v.Fp[n] - v.Fm[n])/(Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(Dx2) - g[2]*v.S[n] +\
                     v.F[n]*(g1p[2] - g1m[2])/(Dx1);
                     F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) + \
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     g[0]*v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1) +\
                     g[0]*v.G[n]*(g2p[1] - g2m[1])/(g[1]*Dx2);
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[2]*Dx2) - v.S[n];
                     //F[n] = g[2]*(v.Fp[n] - v.Fm[n])/(Dx1) + \
                     g[2]*(v.Gp[n] - v.Gm[n])/(Dx2) - g[2]*v.S[n];
                     //F[n] = g[0]*(v.Fp[n] - v.Fm[n])/(Dx1) + \
                     g[0]*(v.Gp[n] - v.Gm[n])/(Dx2) - g[0]*v.S[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[1]*Dx2) - v.S[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[1]*Dx2) - v.S[n]/g[1];
                  }

                  for(n = 0; n < eq; n++)
                  {
                     q1[c2(n,i,j)] = q[c2(n,i,j)] - (Dt)*(F[n]);
                     //q1[c2(n,i,j)] = q[c2(n,i,j)] - (Dt)*(F[n]);
                  }
               }
            }
         break;

         case 2:
            #pragma omp for
            for(i = 4; i <= Nx1-4; i++)
            {
               for(j = 4; j <= Nx2-4; j++)
               {
                  I[0] = i;
                  I[1] = j;
                  r = 'G';

                  RECONST2D(u1,r,&l,I);
                  FLUX2D(&v,&l,I);

                  gauge(g,X1[i],X2[j],0);
                  gauge(g1p,X1p[i],X2[j],0);
                  gauge(g1m,X1m[i],X2[j],0);
                  gauge(g2p,X1[i],X2p[j],0);
                  gauge(g2m,X1[i],X2m[j],0);

                  for(n = 0; n < eq; n++)
                  {
                     //F[n] = g[1]*(g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     g[1]*(g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[2]*v.S[n] +\
                     (g1p[1] - g1m[1])/(Dx1)*v.F[n] +\
                     (g2p[1] - g2m[1])/(Dx2)*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(g[0]*Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(g[0]*Dx2) -\
                     v.S[n] +\
                     (g1p[1] - g1m[1])/(g[1]*Dx1)*v.F[n] +\
                     (g2p[1] - g2m[1])/(g[1]*Dx2)*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     (g1p[1] - g1m[1])/(g[1]*Dx1)*g[0]*v.F[n] +\
                     (g2p[1] - g2m[1])/(g[1]*Dx2)*g[0]*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     (g1p[1] - g1m[1])/(g[1]*Dx1)*g[0]*v.F[n] +\
                     (g2p[1] - g2m[1])/(g[1]*Dx2)*g[0]*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(g[0]*Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(g[0]*Dx2) -\
                     v.S[n] +\
                     (g1p[1] - g1m[1])/(Dx1*g[1])*v.F[n] +\
                     (g2p[1] - g2m[1])/(Dx2*g[1])*v.G[n];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     (g1p[1] - g1m[1])/(Dx1*g[1])*v.F[n]*g[0] +\
                     (g2p[1] - g2m[1])/(Dx2*g[1])*v.G[n]*g[0];
                     //F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) +\
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     (g1p[1] - g1m[1])/(Dx1*g[1])*v.F[n]*g[0] -\
                     (g2p[1] - g2m[1])/(Dx2*g[1])*v.G[n]*g[0];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[2]*Dx2) - v.S[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(Dx2) - g[2]*v.S[n];
                     //F[n] = g[1]*(g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(Dx2) - g[2]*v.S[n] +\
                     g[0]*v.F[n]*(g1p[1] - g1m[1])/(Dx1);
                     //F[n] = g[1]*(g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) + \
                     g[1]*(g2m[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[2]*v.S[n] +\
                     g[0]*v.F[n]*(g1p[1] - g1m[1])/(Dx1) +\
                     g[0]*v.G[n]*(g2p[1] - g2m[1])/(Dx1);
                     F[n] = (g1p[0]*v.Fp[n] - g1m[0]*v.Fm[n])/(Dx1) + \
                     (g2p[0]*v.Gp[n] - g2m[0]*v.Gm[n])/(Dx2) -\
                     g[0]*v.S[n] +\
                     g[0]*v.F[n]*(g1p[1] - g1m[1])/(g[1]*Dx1) +\
                     g[0]*v.G[n]*(g2p[1] - g2m[1])/(g[1]*Dx2);
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[2]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[2]*Dx2) - v.S[n];
                     //F[n] = g[2]*(v.Fp[n] - v.Fm[n])/(Dx1) + \
                     g[2]*(v.Gp[n] - v.Gm[n])/(Dx2) - g[2]*v.S[n];
                     //F[n] = (v.Fp[n] - v.Fm[n])/(g[1]*Dx1) + \
                     (v.Gp[n] - v.Gm[n])/(g[1]*Dx2) - v.S[n]/g[1];
                     //F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
                     (v.Gp[n] - v.Gm[n])/(Dx2) - v.S[n];
                     //F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) + \
                     (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[1]*Dx2) - v.S[n]/g[1];
                  }

                  for(n = 0; n < eq; n++)
                  {
                     q2[c2(n,i,j)] = 0.5*(q1[c2(n,i,j)] + q[c2(n,i,j)] - (Dt)*F[n]);
                     //q2[c2(n,i,j)] = 0.5*(q1[c2(n,i,j)] + q[c2(n,i,j)] - (Dt)*F[n]);
                  }
               }
            }
         break;
      }
   }

   return 0;
}

int RK3D(double *u, double *u1, double *u2, int order)
{
   char r;
   int n, i, j, k;
   int I[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dx3 = dx3;
   double Dt  = dt;
   double L[eq+1], F[eq+1];
   vec_ v;
   lim_ l;

   #pragma omp parallel shared(u,Dt,Dx1,order,r) private(n,L,F,v,l,i,j,k) num_threads(4)
   {
      switch(order)
      {
         case 1:
            #pragma omp for
            for(i = 3; i <= Nx1-3; i++)
            {
               for(j = 3; j <= Nx2-3; j++)
               {
                  for(k = 3; k <= Nx3-3; k++)
                  {
                     //#pragma omp atomic read
                     I[0] = i;
                     I[1] = j;
                     I[2] = k;
                     r = 'G';

                     RECONST3D(u,r,&l,I);
                     FLUX3D(&v,&l,I);

                     for(n = 0; n < eq; n++)
                     {
                        F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
                        (v.Gp[n] - v.Gm[n])/(Dx2) + \
                        (v.Hp[n] - v.Hm[n])/(Dx3) - v.S[n];
                     }

                     MxV(v.A,F,L);

                     for(n = 0; n < eq; n++)
                     {
                        u1[c3(n,i,j,k)] = u[c3(n,i,j,k)] - (dt)*L[n];
                     }
                  }
               }
            }
         break;

         case 2:
            #pragma omp for
            for(i = 3; i <= Nx1-3; i++)
            {
               for(j = 3; j <= Nx2-3; j++)
               {
                  for(k = 3; k <= Nx3-3; k++)
                  {
                     //#pragma omp atomic read
                     I[0] = i;
                     I[1] = j;
                     I[2] = k;
                     r = 'C';

                     RECONST3D(u1,r,&l,I);
                     FLUX3D(&v,&l,I);

                     for(n = 0; n < eq; n++)
                     {
                        F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
                        (v.Gp[n] - v.Gm[n])/(Dx2) + \
                        (v.Hp[n] - v.Hm[n])/(Dx3) - v.S[n];
                     }

                     MxV(v.A,F,L);

                     for(n = 0; n < eq; n++)
                     {
                        u2[c3(n,i,j,k)] = 0.5*(u1[c3(n,i,j,k)] + \
                        u[c3(n,i,j,k)] - (dt)*L[n]);
                     }
                  }
               }
            }
         break;
      }
   }

   return 0;
}


int MxV(double *M, double *V, double *L)
{
   int n, m;
   double res=0.0;

   for(m = 0; m < eq; m++)
   {
      for(n = 0; n < eq; n++)
      {
         res += M[m*eq + n]*V[n];
      }

      L[m] = res;
      res = 0.0;
   }

   return 0;
}

void vel(double *u)
{
   if(*u > 1.0)
   {
      *u = 0.90;
   }
   else if(*u < -1.0)
   {
      *u = -0.90;
   }
}
