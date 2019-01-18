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
      //#pragma omp atomic read
      I[0] = i;
      r = limiter;
 
      RECONST1D(u,r,&l,I);
      FLUX1D(&v,&l,I);
 
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

   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            I[0] = i;
            I[1] = j;
            r = limiter;

            RECONST2D(u,r,&l,I);
            FLUX2D(&v,&l,I);
            
            for(n = 0; n < eq; n++)
            {
               F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
               (v.Gp[n] - v.Gm[n])/(Dx2) - \
               v.S[n];
               //F[n] = roundgen(F[n]);
            }

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
//            if((i == gc && j == Nx2/2-3))
//            {
//               printf("0 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Fp[0],Fm[0],Gp[0],Gm[0],S[0]);
//               printf("1 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Fp[1],Fm[1],Gp[1],Gm[1],S[1]);
//               printf("2 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Fp[2],Fm[2],Gp[2],Gm[2],S[2]);
//               printf("3 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Fp[3],Fm[3],Gp[3],Gm[3],S[3]);
//               printf("0 dF %e dG %e\n",Fp[0]-Fm[0],Gp[0]-Gm[0]);
//               printf("1 dF %e dG %e\n",Fp[1]-Fm[1],Gp[1]-Gm[1]);
//               printf("2 dF %e dG %e\n",Fp[2]-Fm[2],Gp[2]-Gm[2]);
//               printf("3 dF %e dG %e\n",Fp[3]-Fm[3],Gp[3]-Gm[3]);
//               printf("0 F %e \n",F[0]);
//               printf("1 F %e \n",F[1]);
//               printf("2 F %e \n",F[2]);
//               printf("3 F %e \n",F[3]);
//               printf("0 q1 %e \n",q1[c2(0,i,j)]);
//               printf("1 q1 %e \n",q1[c2(1,i,j)]);
//               printf("2 q1 %e \n",q1[c2(2,i,j)]);
//               printf("3 q1 %e \n",q1[c2(3,i,j)]);
//               printf("0 q2 %e \n",q2[c2(0,i,j)]);
//               printf("1 q2 %e \n",q2[c2(1,i,j)]);
//               printf("2 q2 %e \n",q2[c2(2,i,j)]);
//               printf("3 q2 %e \n",q2[c2(3,i,j)]);
//               getchar();                                         
//            }                                                     
//            if((j == gc && i == Nx1/2-3))                         
//            {                                                     
//               printf("0 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Gp[0],Gm[0],Fp[0],Fm[0],S[0]);
//               printf("1 %e %e Fp %e Fm %e Gp %e Gm %e S %e\n",X1[i],X2[j],Gp[1],Gm[1],Fp[1],Fm[1],S[1]);
//               printf("3 %e %e Gp %e Gm %e Fp %e Fm %e S %e\n",X1[i],X2[j],Gp[3],Gm[3],Fp[3],Fm[3],S[3]);
//               printf("2 %e %e Gp %e Gm %e Fp %e Fm %e S %e\n",X1[i],X2[j],Gp[2],Gm[2],Fp[2],Fm[2],S[2]);
//               printf("0 dF %e dG %e\n",Fp[0]-Fm[0],Gp[0]-Gm[0]);
//               printf("1 dF %e dG %e\n",Fp[1]-Fm[1],Gp[1]-Gm[1]);
//               printf("3 dF %e dG %e\n",Fp[3]-Fm[3],Gp[3]-Gm[3]);
//               printf("2 dF %e dG %e\n",Fp[2]-Fm[2],Gp[2]-Gm[2]);
//               printf("0 F %e \n",F[0]);
//               printf("1 F %e \n",F[1]);
//               printf("2 F %e \n",F[2]);
//               printf("3 F %e \n",F[3]);
//               printf("0 q1 %e \n",q1[c2(0,i,j)]);
//               printf("1 q1 %e \n",q1[c2(1,i,j)]);
//               printf("2 q1 %e \n",q1[c2(2,i,j)]);
//               printf("3 q1 %e \n",q1[c2(3,i,j)]);
//               printf("0 q2 %e \n",q2[c2(0,i,j)]);
//               printf("1 q2 %e \n",q2[c2(1,i,j)]);
//               printf("2 q2 %e \n",q2[c2(2,i,j)]);
//               printf("3 q2 %e \n",q2[c2(3,i,j)]);
//               getchar();
//            }
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
               r = limiter;

               RECONST3D(u,r,&l,I);
               FLUX3D(&v,&l,I);

               for(n = 0; n < eq; n++)
               {
                  F[n] = (v.Fp[n] - v.Fm[n])/(Dx1) + \
                  (v.Gp[n] - v.Gm[n])/(Dx2) + \
                  (v.Hp[n] - v.Gm[n])/(Dx3) - \
                  v.S[n];
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
