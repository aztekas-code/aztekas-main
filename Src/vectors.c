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

int AMATRIX1D(double *u, vec_ *v, int *I)
{
   int m, n;
   int i;
   double geoS[eq + 1], extS[eq + 1];

   x1 = grid.X1[I[0]];
#if COORDINATES == 0 || COORDINATES == 1
   x2 = 0;
#elif COORDINATES == 2
   x2 = M_PI_2;
#endif
   x3 = 0;

   Source_Terms(geoS,u);
   User_Source_Terms(extS,u);

#if integration == 1
   funct_A(v->A,u);
#endif

   for(n = 0; n < eq; n++)
   {
      v->S[n] = geoS[n] + extS[n];
   }

   return 0;
}

int AMATRIX2D(double *u, vec_ *v, int *I)
{
   int m, n;
   double geoS[eq + 1], extS[eq + 1];

   x1  = grid.X1[I[0]];
   x2  = grid.X2[I[1]];
   x3  = 0;

#if polar == 1
   x2 = M_PI_2;
#endif
   
   Source_Terms(geoS,u);
   User_Source_Terms(extS,u);

#if integration == 1
   funct_A(v->A,u);
#endif

   for(n = 0; n < eq; n++)
   {
      v->S[n] = geoS[n] + extS[n];
   }

   return 0;
}

int AMATRIX3D(double *u, vec_ *v, int *I)
{
   int m, n;

   x1  = grid.X1[I[0]];
   x2  = grid.X2[I[1]];
   x3  = grid.X3[I[2]];

   return 0;
}

///////////////////////////////////////////////////////////////////////////

int VECTOR(int pm, char flux, lim_ *l, flx_ *f, int *I)
{
   int n;
   double *u, lr, ll;
   double up[eq + 1];
   double um[eq + 1];
   double dp[3];
   double dm[3];
   double dup[eq + 1];
   double dum[eq + 1];
   x1 = 0.0;
   x2 = 0.0;
   x3 = 0.0;

#if DIM == 1
   x1 = grid.X1[I[0]];
   x2 = M_PI_2;
#elif DIM == 2
   x1 = grid.X1[I[0]];
   x2 = grid.X2[I[1]];
#elif DIM == 3 
   x1 = grid.X1[I[0]];
   x2 = grid.X2[I[1]];
   x3 = grid.X3[I[2]];
#endif

   if(pm == 1)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1p;
            x1 = grid.X1p[I[0]];
         break;

         case 'g':
            u = l->ux2p;
            x2 = grid.X2p[I[1]];
         break;

         case 'h':
            u = l->ux3p;
            x3 = grid.X3p[I[2]];
         break;
      }
   }
   if(pm == 0)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1m;
            x1 = grid.X1m[I[0]];
         break;

         case 'g':
            u = l->ux2m;
            x2 = grid.X2m[I[1]];
         break;

         case 'h':
            u = l->ux3m;
            x3 = grid.X3m[I[2]];
         break;
      }
   }

   for(n = 0; n < eq; n++)
   {
      f->up[n] = u[1*eq + n];
      f->um[n] = u[0*eq + n];
   }

   Prim2Cons(f->qp,f->up);
   Prim2Cons(f->qm,f->um);

   switch(flux)
   {
      case 'f':
         Prim2FluxF(f->fp,f->up);
         funct_Dm(dp,f->up);

         Prim2FluxF(f->fm,f->um);
         funct_Dm(dm,f->um);
      break;

      case 'g':
         Prim2FluxG(f->fp,f->up);
         funct_Dn(dp,f->up);

         Prim2FluxG(f->fm,f->um);
         funct_Dn(dm,f->um);
      break;

      case 'h':
         Prim2FluxH(f->fp,f->up);
         funct_Do(dp,f->up);

         Prim2FluxH(f->fm,f->um);
         funct_Do(dm,f->um);
      break;
   }

   lr = MAX(dp[0],dp[1]);
   lr = MAX(lr,dp[2]);
   lr = MAX(0.0,lr);
   ll = MAX(dm[0],dm[1]);
   ll = MAX(ll,dm[2]);
   ll = MAX(0.0,ll);

   f->lp = MAX(lr,ll);

   lr = MIN(dp[0],dp[1]);
   lr = MIN(lr,dp[2]);
   lr = MIN(0.0,lr);
   ll = MIN(dm[0],dm[1]);
   ll = MIN(ll,dm[2]);
   ll = MIN(0.0,ll);

   f->lm = MIN(lr,ll);

   return 0;
}

///////////////////////////////////////////////////////////////////////////
