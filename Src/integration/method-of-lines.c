/**
 * @file /integration/runge-kutta.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Principal loop in the integration.
 */

#include"main.h"

#if DIM == 1

void Method_of_Lines(int order)
{
   int cell[3];
   double Dx1 = dx1;
   double Dt  = dt;
   rhs_ rhs;
   rk_ rk;

#ifdef _OPENMP
   #pragma omp parallel default(none) \
               shared(U,Q,Q0,Q1,Q2,grid,Nx1,Dt,order,U1p,U1m)
   #pragma omp for private(cell,rk,rhs,Dx1)
#endif
   for(int i = gc; i <= Nx1-gc; i++)
   {
      cell[0] = i;
 
      Dx1 = grid.X1p[i] - grid.X1m[i];

      Numerical_Flux_F(rhs.Fp,PLUS,cell);
      Numerical_Flux_F(rhs.Fm,MINUS,cell);

      Prim2Sources(rhs.S,cell);

      for(int n = 0; n < eq; n++)
      {
         rhs.F[n] = - (S1p(i)*rhs.Fp[n] - S1m(i)*rhs.Fm[n])/(Dx1) + \
                      rhs.S[n];
      }
 
#if INTEGRATION == PVRS //PVRS
      for(int n = 0; n < eq; n++)
      {
         rhs.L[n] = rhs.F[n];
      }

      MxV(rhs.A,rhs.L,rhs.F);
#endif

      switch(order)
      {
         case 1:
            for(int n = 0; n < eq; n++)
            {
               rk.u0 = Q0(n,i);
               rk.f  = rhs.F[n];
               rk.h  = Dt;

               Runge_Kutta(&rk,1);

               Q0(n,i) = rk.u0;
               Q1(n,i) = rk.u1;
               Q(n,i)  = rk.u1;
            }
         break;
 
         case 2:
            for(int n = 0; n < eq; n++)
            {
               rk.u0 = Q0(n,i);
               rk.u1 = Q1(n,i);
               rk.f  = rhs.F[n];
               rk.h  = Dt;

               Runge_Kutta(&rk,2);

               Q(n,i) = rk.u2;
            }
         break;
      }
   }
}

#elif DIM == 2 || DIM == 4

void Method_of_Lines(int order)
{
   int cell[3];
   double Dx1 = dx1;
   double Dx2 = dx2;
   double Dt  = dt;
   rhs_ rhs;

#ifdef _OPENMP
   #pragma omp parallel default(none) \
               shared(U,Q,Q0,Q1,Q2,grid,Nx1,Nx2,Dt,order,U1p,U1m)
   #pragma omp for private(cell,rhs,Dx1,Dx2) collapse(2)
#endif
   for(int j = gc; j <= Nx2-gc; j++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         cell[0] = i;
         cell[1] = j;

         Dx1 = grid.X1p[i] - grid.X1m[i];
         Dx2 = grid.X2p[j] - grid.X2m[j];

         Numerical_Flux_F(rhs.Fp,PLUS,cell);
         Numerical_Flux_F(rhs.Fm,MINUS,cell);

         Numerical_Flux_G(rhs.Gp,PLUS,cell);
         Numerical_Flux_G(rhs.Gm,MINUS,cell);

         Prim2Sources(rhs.S,cell);

         for(int n = 0; n < eq; n++)
         {
            rhs.F[n] = (S1p(i,j)*rhs.Fp[n] - S1m(i,j)*rhs.Fm[n])/(Dx1) + \
                       (S2p(i,j)*rhs.Gp[n] - S2m(i,j)*rhs.Gm[n])/(Dx2) - \
                       rhs.S[n];
         }

#if INTEGRATION == PVRS //PVRS
         for(int n = 0; n < eq; n++)
         {
            rhs.L[n] = rhs.F[n];
         }

         MxV(rhs.A,rhs.L,rhs.F);
#endif

         switch(order)
         {
            case 1:
               for(int n = 0; n < eq; n++)
               {
                  Q1(n,i,j) = Q0(n,i,j) - (Dt)*(rhs.F[n]);
                  Q(n,i,j)  = Q1(n,i,j);
               }
            break;

            case 2:
               for(int n = 0; n < eq; n++)
               {
                  Q2(n,i,j) = 0.5*(Q1(n,i,j) + Q0(n,i,j) - (Dt)*(rhs.F[n]));
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
   rhs_ v;
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
               F[n] = (S1p(i,j,k)*rhs.Fp[n] - S1m(i,j,k)*rhs.Fm[n])/(Dx1) + \
                      (S2p(i,j,k)*rhs.Gp[n] - S2m(i,j,k)*rhs.Gm[n])/(Dx2) - \
                      (S3p(i,j,k)*rhs.Hp[n] - S3m(i,j,k)*rhs.Hm[n])/(Dx3) - \
                      rhs.S[n];
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
