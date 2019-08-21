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
#include"main.h"

#if DIM == 1

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
 
      Reconst1D(u,&l,I);
      Flux1D(&v,&l,I);
 
      for(n = 0; n < eq; n++)
      {
         F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) - \
         v.S[n];
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
               q1(n,i) = q(n,i) - (Dt)*(F[n]);
            }
         break;

         case 2:
            for(n = 0; n < eq; n++)
            {
               q2(n,i) = 0.5*(q1(n,i) + q(n,i) - (Dt)*F[n]);
            }
         break;
      }
   }

   return 0;
}

#elif DIM == 2 || DIM == 4

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

         Dx1 = grid.X1p[i] - grid.X1m[i];

         Reconst2D(u,&l,I);
         Flux2D(&v,&l,I);

         for(n = 0; n < eq; n++)
         {
            F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
            (v.Gp[n] - v.Gm[n])/(grid.X1[i]*Dx2) - \
            v.S[n];
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
                  q1(n,i,j) = q(n,i,j) - (Dt)*(F[n]);
               }
            break;

            case 2:
               for(n = 0; n < eq; n++)
               {
                  q2(n,i,j) = 0.5*(q1(n,i,j) + q(n,i,j) - (Dt)*F[n]);
               }
            break;
         }
      }
   }

   return 0;
}

#elif DIM == 3

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

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            I[0] = i;
            I[1] = j;
            I[2] = k;

            Reconst3D(u,&l,I);
            Flux3D(&v,&l,I);

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
                     q1(n,i,j,k) = q(n,i,j,k) - (Dt)*L[n];
                  }
               break;

               case 2:
                  for(n = 0; n < eq; n++)
                  {
                     q2(n,i,j,k) = 0.5*(q1(n,i,j,k) + \
                     q(n,i,j,k) - (Dt)*L[n]);
                  }
               break;
            }
         }
      }
   }

   return 0;
}

#endif
