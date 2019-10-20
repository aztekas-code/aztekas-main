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

int Reconst1D(double *u, lim_ *l, int *I)
{
   int n, i;
   int reconst;
   double dup2, dup1, dum1, dum2;

   i = I[0];
   reconst = RECONST;
#if x1min_exc == TRUE
   if (i <= 2*gc && RECONST == WENO5)
   {
      reconst = MC;
   }
   if (i <= 2*gc && RECONST != WENO5)
   {
      reconst = GODUNOV;
   }
#endif

   for(n = 0; n < eq; n++)
   {
      l->ux[n] = u(n,i);

      if(reconst == WENO5)
      {
         // U1m(i-1)
         l->ux1m[0*eq + n] = Weno5(u(n,i-3),\
                                   u(n,i-2),\
                                   u(n,i-1),\
                                   u(n, i ),\
                                   u(n,i+1));
         // U1p(i)
         l->ux1m[1*eq + n] = Weno5(u(n,i+2),\
                                   u(n,i+1),\
                                   u(n, i ),\
                                   u(n,i-1),\
                                   u(n,i-2));

         // U1m(i)
         l->ux1p[0*eq + n] = Weno5(u(n,i-2),\
                                   u(n,i-1),\
                                   u(n, i ),\
                                   u(n,i+1),\
                                   u(n,i+2));
         // U1p(i+1)
         l->ux1p[1*eq + n] = Weno5(u(n,i+3),\
                                   u(n,i+2),\
                                   u(n,i+1),\
                                   u(n, i ),\
                                   u(n,i-1));
      }
      else
      {
         dum2 = u(n,i-1) - u(n,i-2);
         dum1 = u(n, i ) - u(n,i-1);
         dup1 = u(n,i+1) - u(n, i );
         dup2 = u(n,i+2) - u(n,i+1);

         l->sx1[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx1[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx1[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux1m[0*eq + n] = u(n,i-1) + 0.5*l->sx1[0*eq + n];
         l->ux1m[1*eq + n] = u(n, i ) - 0.5*l->sx1[1*eq + n];

         l->ux1p[0*eq + n] = u(n, i ) + 0.5*l->sx1[1*eq + n];
         l->ux1p[1*eq + n] = u(n,i+1) - 0.5*l->sx1[2*eq + n];
      }
   }

   return 0;
}

#elif DIM == 2 || DIM == 4

int Reconst2D(double *u, lim_ *l, int *I)
{
   int n, i, j;
   int reconst;
   double dup2, dup1, dum1, dum2;

   i = I[0];
   j = I[1];
   reconst = RECONST;

#if x1min_exc == TRUE
   if (i <= 2*gc && RECONST == WENO5)
   {
      reconst = MC;
   }
   if (i <= 2*gc && RECONST != WENO5)
   {
      reconst = GODUNOV;
   }
#endif

   for(n = 0; n < eq; n++)
   {
      l->ux[n] = u(n,i,j);

      if(reconst == WENO5)
      {
         l->ux1m[0*eq + n] = Weno5(u(n,i-3,j),\
                                   u(n,i-2,j),\
                                   u(n,i-1,j),\
                                   u(n, i ,j),\
                                   u(n,i+1,j));
         l->ux1m[1*eq + n] = Weno5(u(n,i+2,j),\
                                   u(n,i+1,j),\
                                   u(n, i ,j),\
                                   u(n,i-1,j),\
                                   u(n,i-2,j));

         l->ux1p[0*eq + n] = Weno5(u(n,i-2,j),\
                                   u(n,i-1,j),\
                                   u(n, i ,j),\
                                   u(n,i+1,j),\
                                   u(n,i+2,j));
         l->ux1p[1*eq + n] = Weno5(u(n,i+3,j),\
                                   u(n,i+2,j),\
                                   u(n,i+1,j),\
                                   u(n, i ,j),\
                                   u(n,i-1,j));

         l->ux2m[0*eq + n] = Weno5(u(n,i,j-3),\
                                   u(n,i,j-2),\
                                   u(n,i,j-1),\
                                   u(n,i, j ),\
                                   u(n,i,j+1));
         l->ux2m[1*eq + n] = Weno5(u(n,i,j+2),\
                                   u(n,i,j+1),\
                                   u(n,i, j ),\
                                   u(n,i,j-1),\
                                   u(n,i,j-2));

         l->ux2p[0*eq + n] = Weno5(u(n,i,j-2),\
                                   u(n,i,j-1),\
                                   u(n,i, j ),\
                                   u(n,i,j+1),\
                                   u(n,i,j+2));
         l->ux2p[1*eq + n] = Weno5(u(n,i,j+3),
                                   u(n,i,j+2),\
                                   u(n,i,j+1),\
                                   u(n,i, j ),\
                                   u(n,i,j-1));
      }
      else
      {
         dum2 = u(n,i-1,j) - u(n,i-2,j);
         dum1 = u(n, i ,j) - u(n,i-1,j);
         dup1 = u(n,i+1,j) - u(n, i ,j);
         dup2 = u(n,i+2,j) - u(n,i+1,j);

         l->sx1[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx1[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx1[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux1m[0*eq + n] = u(n,i-1,j) + 0.5*l->sx1[0*eq + n];
         l->ux1m[1*eq + n] = u(n, i ,j) - 0.5*l->sx1[1*eq + n];

         l->ux1p[0*eq + n] = u(n, i ,j) + 0.5*l->sx1[1*eq + n];
         l->ux1p[1*eq + n] = u(n,i+1,j) - 0.5*l->sx1[2*eq + n];

         dum2 = u(n,i,j-1) - u(n,i,j-2);
         dum1 = u(n,i, j ) - u(n,i,j-1);
         dup1 = u(n,i,j+1) - u(n,i, j );
         dup2 = u(n,i,j+2) - u(n,i,j+1);

         l->sx2[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx2[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx2[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux2m[0*eq + n] = u(n,i,j-1) + 0.5*l->sx2[0*eq + n];
         l->ux2m[1*eq + n] = u(n,i,j  ) - 0.5*l->sx2[1*eq + n];

         l->ux2p[0*eq + n] = u(n,i,j  ) + 0.5*l->sx2[1*eq + n];
         l->ux2p[1*eq + n] = u(n,i,j+1) - 0.5*l->sx2[2*eq + n];
      }
   }

   return 0;
}

#elif DIM == 3

int Reconst3D(double *u, lim_ *l, int *I)
{
   int n, i, j, k;
   int reconst;
   double dup2, dup1, dum1, dum2;

   i = I[0];
   j = I[1];
   k = I[2];
   reconst = RECONST;

   for(n = 0; n < eq; n++)
   {
      l->ux[n] = u(n,i,j,k);

      if(reconst == WENO5)
      {
         l->ux1m[0*eq + n] = Weno5(u(n,i-3,j,k),\
                                   u(n,i-2,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i+1,j,k));
         l->ux1m[1*eq + n] = Weno5(u(n,i+2,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n,i-2,j,k));

         l->ux1p[0*eq + n] = Weno5(u(n,i-2,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n,i+2,j,k));
         l->ux1p[1*eq + n] = Weno5(u(n,i+3,j,k),\
                                   u(n,i+2,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i-1,j,k));

         l->ux2m[0*eq + n] = Weno5(u(n,i,j-3,k),\
                                   u(n,i,j-2,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j+1,k));
         l->ux2m[1*eq + n] = Weno5(u(n,i,j+2,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i,j-2,k));

         l->ux2p[0*eq + n] = Weno5(u(n,i,j-2,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i,j+2,k));
         l->ux2p[1*eq + n] = Weno5(u(n,i,j+3,k),\
                                   u(n,i,j+2,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j-1,k));

         l->ux3m[0*eq + n] = Weno5(u(n,i,j,k-3),\
                                   u(n,i,j,k-2),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k+1));
         l->ux3m[1*eq + n] = Weno5(u(n,i,j,k+2),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j,k-2));

         l->ux3p[0*eq + n] = Weno5(u(n,i,j,k-2),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j,k+2));
         l->ux3p[1*eq + n] = Weno5(u(n,i,j,k+3),\
                                   u(n,i,j,k+2),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k-1));
      }
      else
      {
         dum2 = u(n,i-1,j,k) - u(n,i-2,j,k);
         dum1 = u(n,i  ,j,k) - u(n,i-1,j,k);
         dup1 = u(n,i+1,j,k) - u(n,i  ,j,k);
         dup2 = u(n,i+2,j,k) - u(n,i+1,j,k);

         l->sx1[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx1[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx1[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux1m[0*eq + n] = u(n,i-1,j,k) + 0.5*l->sx1[0*eq + n];
         l->ux1m[1*eq + n] = u(n,i  ,j,k) - 0.5*l->sx1[1*eq + n];

         l->ux1p[0*eq + n] = u(n,i  ,j,k) + 0.5*l->sx1[1*eq + n];
         l->ux1p[1*eq + n] = u(n,i+1,j,k) - 0.5*l->sx1[2*eq + n];

         dum2 = u(n,i,j-1,k) - u(n,i,j-2,k);
         dum1 = u(n,i,j  ,k) - u(n,i,j-1,k);
         dup1 = u(n,i,j+1,k) - u(n,i,j  ,k);
         dup2 = u(n,i,j+2,k) - u(n,i,j+1,k);

         l->sx2[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx2[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx2[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux2m[0*eq + n] = u(n,i,j-1,k) + 0.5*l->sx2[0*eq + n];
         l->ux2m[1*eq + n] = u(n,i,j  ,k) - 0.5*l->sx2[1*eq + n];

         l->ux2p[0*eq + n] = u(n,i,j  ,k) + 0.5*l->sx2[1*eq + n];
         l->ux2p[1*eq + n] = u(n,i,j+1,k) - 0.5*l->sx2[2*eq + n];

         dum2 = u(n,i,j,k-1) - u(n,i,j,k-2);
         dum1 = u(n,i,j,k  ) - u(n,i,j,k-1);
         dup1 = u(n,i,j,k+1) - u(n,i,j,k  );
         dup2 = u(n,i,j,k+2) - u(n,i,j,k+1);

         l->sx3[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx3[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx3[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux3m[0*eq + n] = u(n,i,j,k-1) + 0.5*l->sx3[0*eq + n];
         l->ux3m[1*eq + n] = u(n,i,j,k  ) - 0.5*l->sx3[1*eq + n];

         l->ux3p[0*eq + n] = u(n,i,j,k  ) + 0.5*l->sx3[1*eq + n];
         l->ux3p[1*eq + n] = u(n,i,j,k+1) - 0.5*l->sx3[2*eq + n];
      }
   }

   return 0;
}

#endif

double Limiter(double A, double B, int reconst)
{
   double sig;

   switch(reconst)
   {
      case GODUNOV:
         sig = Godunov(A,B);
      break;

      case MINMOD:
         sig = Minmod(A,B);
      break;

      case MC:
         sig = Mc(A,B);
      break;

      case SUPERBEE:
         sig = Superbee(A,B);
      break;
   }

   return sig;
}

double Maxmod(double A, double B)
{
   double sig;

   if(A*B < 0)
   {
      sig = 0;
   }
   else if(A*B >= 0)
   {
      if(fabs(A) > fabs(B))
      {
         sig = A;
      }
      else if(fabs(A) < fabs(B))
      {
         sig = B;
      }
   }

   return sig;
}

double Minmod(double A, double B)
{
   double sig;

   if(A*B <= 0)
   {
      sig = 0;
   }
   else if(A*B > 0)
   {
      if(fabs(A) < fabs(B))
      {
         sig = A;
      }
      else if(fabs(A) > fabs(B))
      {
         sig = B;
      }
   }

   return sig;
}

double Mc(double A, double B)
{
   double sig;
   double C;

   C = (A + B) / 2.0e+00;

   if(A*B <= 0)
   {
      sig = 0.0e+00;
   }
   else if(A*B > 0)
   {
      if(fabs(A) < fabs(B) && 2.0e+00*fabs(A) <= fabs(C))
      {
         sig = 2.0e+00*A;
      }
      else if(fabs(A) >= fabs(B) && 2.0e+00*fabs(B) <= fabs(C))
      {
         sig = 2.0e+00*B;
      }
      else if(fabs(C) <= 2.0e+00*fabs(A) && fabs(C) <= 2.0e+00*fabs(B))
      {
         sig= C;
      }
   }

   return sig;
}

double Superbee(double A, double B)
{
   double sig1;
   double sig2;
   double sig;

   sig1 = Minmod(2*A,B);
   sig2 = Minmod(A,2*B);

   sig = Maxmod(sig1,sig2);
   return sig;
}

double Godunov(double A, double B)
{
   return 0.0;
}

double Weno5(double v1, double v2, double v3, double v4, double v5)
{
   double c13d12 = 1.3e+01 / 1.2e+01;
   double c1d4   = 1.0e+00 / 4.0e+00;
   double c11d6  = 1.1e+01 / 6.0e+00;
   double c2d6  = 2.0e+00 / 6.0e+00;
   double c7d6  = 7.0e+00 / 6.0e+00;
   double c5d6  = 5.0e+00 / 6.0e+00;
   double c1d6  = 1.0e+00 / 6.0e+00;
   double c1d3  = 1.0e+00 / 3.0e+00;

   double fs1, fs2, fs3;
   double bs1, bs2, bs3;

   double ep = 1.0e-06;

   double P1, P2, P3;
   double s1, s2, s3;
   double a1, a2, a3;

   double asum;
   double g1, g2, g3;
   double d1, d2, d3;

   double om1, om2, om3;
   double R;

   d1 = 1.0e-01;
   d2 = 6.0e-01;
   d3 = 3.0e-01;

   P1 =  c1d3*v1 - c7d6*v2 + c11d6*v3;
   P2 = -c1d6*v2 + c5d6*v3 +  c1d3*v4;
   P3 =  c1d3*v3 + c5d6*v4 -  c1d6*v5;

   fs1 = (v1 - 2.0e+00 * v2 + v3);
   fs2 = (v2 - 2.0e+00 * v3 + v4);
   fs3 = (v3 - 2.0e+00 * v4 + v5);

   bs1 = (v1 - 4.0e+00*v2 + 3.0e+00*v3);
   bs2 = (v2 - v4);
   bs3 = (3.0e+00*v3 - 4.0e+00*v4 + v5);

   s1 = c13d12*fs1*fs1 + c1d4*bs1*bs1;
   s2 = c13d12*fs2*fs2 + c1d4*bs2*bs2;
   s3 = c13d12*fs3*fs3 + c1d4*bs3*bs3;

   a1 = d1 * (1.0e+00 / pow(s1 + ep,2.0e+00));
   a2 = d2 * (1.0e+00 / pow(s2 + ep,2.0e+00));
   a3 = d3 * (1.0e+00 / pow(s3 + ep,2.0e+00));

   asum = a1 + a2 + a3;

   om1 = a1 / asum;
   om2 = a2 / asum;
   om3 = a3 / asum;

   R = om1*P1 + om2*P2 + om3*P3;

   return R;
}
