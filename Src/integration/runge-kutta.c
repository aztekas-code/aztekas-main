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

void Runge_Kutta(int order)
{
   int cell[3];
   double Dx1 = dx1;
   double Dt  = dt;
   rhs_ vec;

   #pragma omp parallel default(none) shared(U,Q,Q0,Q1,Q2,grid,Nx1,Dt,order,U1p,U1m) if (OMP_NUM > 1)
   #pragma omp for private(cell,vec,Dx1) 
   for(int i = gc; i <= Nx1-gc; i++)
   {
      cell[0] = i;
 
      Dx1 = grid.X1p[i] - grid.X1m[i];

      Numerical_Flux_F(vec.Fp,PLUS,cell);
      Numerical_Flux_F(vec.Fm,MINUS,cell);

      Prim2Sources(vec.S,cell);

      for(int n = 0; n < eq; n++)
      {
         vec.F[n] = (S1p(i)*vec.Fp[n] - S1m(i)*vec.Fm[n])/(Dx1) - \
                    vec.S[n];
      }
 
#if INTEGRATION == PVRS //PVRS
     for(int n = 0; n < eq; n++)
     {
        vec.L[n] = vec.F[n];
     }

     MxV(vec.A,vec.L,vec.F);
#endif

      switch(order)
      {
         case 1:
            for(int n = 0; n < eq; n++)
            {
               Q1(n,i) = Q0(n,i) - (Dt)*(vec.F[n]);
               Q(n,i)  = Q1(n,i);
            }
         break;

         case 2:
            for(int n = 0; n < eq; n++)
            {
               Q2(n,i) = 0.5*(Q1(n,i) + Q0(n,i) - (Dt)*(vec.F[n]));
               Q(n,i)  = Q2(n,i);
            }
         break;
      }
   }
}

#elif DIM == 2 || DIM == 4

void Runge_Kutta(int order)
{
   int cell[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dt  = dt;
   rhs_ vec;

   #pragma omp parallel default(none) shared(U,Q,Q0,Q1,Q2,grid,Nx1,Nx2,Dt,order,U1p,U1m) if (OMP_NUM > 1)
   #pragma omp for private(cell,vec,Dx1,Dx2) collapse(2)
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         cell[0] = i;
         cell[1] = j;
    
         Dx1 = grid.X1p[i] - grid.X1m[i];
         Dx2 = grid.X2p[j] - grid.X2m[j];

         Numerical_Flux_F(vec.Fp,PLUS,cell);
         Numerical_Flux_F(vec.Fm,MINUS,cell);

         Numerical_Flux_G(vec.Gp,PLUS,cell);
         Numerical_Flux_G(vec.Gm,MINUS,cell);

         Prim2Sources(vec.S,cell);

         for(int n = 0; n < eq; n++)
         {
            vec.F[n] = (S1p(i,j)*vec.Fp[n] - S1m(i,j)*vec.Fm[n])/(Dx1) + \
                       (S2p(i,j)*vec.Gp[n] - S2m(i,j)*vec.Gm[n])/(Dx2) - \
                       vec.S[n];

         }
    
#if INTEGRATION == PVRS //PVRS
         for(int n = 0; n < eq; n++)
         {
            vec.L[n] = vec.F[n];
         }

         MxV(vec.A,vec.L,vec.F);
#endif

         switch(order)
         {
            case 1:
               for(int n = 0; n < eq; n++)
               {
                  Q1(n,i,j) = Q0(n,i,j) - (Dt)*(vec.F[n]);
                  Q(n,i,j)  = Q1(n,i,j);
               }
            break;

            case 2:
               for(int n = 0; n < eq; n++)
               {
                  Q2(n,i,j) = 0.5*(Q1(n,i,j) + Q0(n,i,j) - (Dt)*(vec.F[n]));
                  Q(n,i,j)  = Q2(n,i,j);
               }
            break;
         }
      }
   } 
}

#elif DIM == 3

int RK3D(double *u, double *q, double *q1, double *q2, int order)
{
   char r;
   int n, i, j, k;
   int cell[3];
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
            cell[0] = i;
            cell[1] = j;
            cell[2] = k;

            Dx1 = grid.X1p[i] - grid.X1m[i];
            Dx2 = grid.X2p[j] - grid.X2m[j];
            Dx3 = grid.X2p[k] - grid.X2m[k];

            Reconst3D(u,&l,cell);
            Flux3D(&v,&l,cell);

            for(n = 0; n < eq; n++)
            {
               F[n] = (S1p(i,j,k)*vec.Fp[n] - S1m(i,j,k)*vec.Fm[n])/(Dx1) + \
                      (S2p(i,j,k)*vec.Gp[n] - S2m(i,j,k)*vec.Gm[n])/(Dx2) - \
                      (S3p(i,j,k)*vec.Hp[n] - S3m(i,j,k)*vec.Hm[n])/(Dx3) - \
                      vec.S[n];
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
