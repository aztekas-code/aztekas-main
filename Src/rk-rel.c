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
#include"main.h"
#include"vector.h"
#include"param.h"

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
   vec_ v;
   lim_ l;

   for(i = gc; i <= Nx1-gc; i++)
   {
      I[0] = i;
      r = limiter;
 
      RECONST1D(u,r,&l,I);
      FLUX1D(&v,&l,I);
 
      GAUGE(g,X1[i],0,0);
      GAUGE(g1p,X1p[i],0,0);
      GAUGE(g1m,X1m[i],0,0);
 
      for(n = 0; n < eq; n++)
      {
         F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) - \
         g[2]*v.S[n]/g[1];
      }

#if integration == 1 //PVRS
      for(n = 0; n < eq; n++)
      {
         L[n] = F[n];
      }

      MxV(v.A,L,F);
#endif

      switch(order)
      {
         case 1:
            for(n = 0; n < eq; n++)
            {
               q1[c1(n,i)] = q[c1(n,i)] - (Dt)*(F[n]);
            }
         break;

         case 2:
            for(n = 0; n < eq; n++)
            {
               q2[c1(n,i)] = 0.5*(q1[c1(n,i)] + q[c1(n,i)] - (Dt)*F[n]);
            }

         break;
      }
   }

   return 0;
}

int RK2D(double *u, double *q, double *q1, double *q2, int order)
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

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         I[0] = i;
         I[1] = j;
         r = limiter;

         RECONST2D(u,r,&l,I);
         FLUX2D(&v,&l,I);

         GAUGE(g,X1[i],X2[j],0);
         GAUGE(g1p,X1p[i],X2[j],0);
         GAUGE(g1m,X1m[i],X2[j],0);
         GAUGE(g2p,X1[i],X2p[j],0);
         GAUGE(g2m,X1[i],X2m[j],0);

         for(n = 0; n < eq; n++)
         {
            F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) +\
            (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[1]*Dx2) -\
            g[2]*v.S[n]/g[1];
         }

#if integration == 1 //PVRS
         for(n = 0; n < eq; n++)
         {
            L[n] = F[n];
         }

         MxV(v.A,L,F);
#endif

         switch(order)
         {
            case 1:
               for(n = 0; n < eq; n++)
               {
                  q1[c2(n,i,j)] = q[c2(n,i,j)] - (Dt)*(F[n]);
               }
            break;

            case 2:
               for(n = 0; n < eq; n++)
               {
                  q2[c2(n,i,j)] = 0.5*(q1[c2(n,i,j)] + q[c2(n,i,j)] - (Dt)*F[n]);
               }
            break;
         }
      }
   }

   return 0;
}

int RK3D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i, j, k;
   int I[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dx3 = dx3;
   double Dt  = dt;
   double g[6], g1p[6], g1m[6];
   double g2p[6], g2m[6];
   double g3p[6], g3m[6];
   double L[eq+1], F[eq+1];
   vec_ v;
   lim_ l;

   #pragma omp parallel shared(u,Dt,Dx1,order,r) private(n,L,F,v,l,i,j,k) num_threads(4)
   {
      #pragma omp for
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            for(k = gc; k <= Nx3-gc; k++)
            {
               I[0] = i;
               I[1] = j;
               I[2] = k;
               r = 'G';

               GAUGE(g,X1[i],X2[j],X3[k]);
               GAUGE(g1p,X1p[i],X2[j],X3[k]);
               GAUGE(g1m,X1m[i],X2[j],X3[k]);
               GAUGE(g2p,X1[i],X2p[j],X3[k]);
               GAUGE(g2m,X1[i],X2m[j],X3[k]);
               GAUGE(g3p,X1[i],X2[j],X3p[k]);
               GAUGE(g3m,X1[i],X2[j],X3m[k]);
   
               RECONST3D(u,r,&l,I);
               FLUX3D(&v,&l,I);

               for(n = 0; n < eq; n++)
               {
                  F[n] = (g1p[2]*v.Fp[n] - g1m[2]*v.Fm[n])/(g[1]*Dx1) + \
                  (g2p[2]*v.Gp[n] - g2m[2]*v.Gm[n])/(g[1]*Dx2) + \
                  (g3p[2]*v.Hp[n] - g3m[2]*v.Gm[n])/(g[1]*Dx3) - \
                  g[2]*v.S[n]/g[1];
               }

               switch(order)
               {
                  case 1:
                     for(n = 0; n < eq; n++)
                     {
                        q1[c3(n,i,j,k)] = q[c3(n,i,j,k)] - (Dt)*L[n];
                     }
                  break;

                  case 2:
                     for(n = 0; n < eq; n++)
                     {
                        q2[c3(n,i,j,k)] = 0.5*(q1[c3(n,i,j,k)] + \
                        q[c3(n,i,j,k)] - (Dt)*L[n]);
                     }
                  break;
               }
            }
         }
      }
   }

   return 0;
}
