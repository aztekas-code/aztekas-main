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

void Sources(double *u, vec_ *v, int *I)
{
   int n;
   double x[4];
   double default_S[eq+1], user_S[eq+1];
   grid_ local_grid;

   local_grid.x[0] = grid.time;

#if DIM == 1

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = 0.0;
   local_grid.x[3] = 0.0;
   #if COORDINATES == SPHERICAL
   local_grid.x[2] = M_PI_2;
   #endif

#elif DIM == 2 || DIM == 4

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = grid.X2[I[1]];
   local_grid.x[3] = 0.0;
   #if POLAR == TRUE
   local_grid.x[2] = M_PI_2;
   #endif

#elif DIM == 3

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = grid.X2[I[1]];
   local_grid.x[3] = grid.X3[I[2]];

#endif

   Source_Terms(default_S,u,local_grid);
   User_Source_Terms(user_S,u,local_grid);

   for(n = 0; n < eq; n++)
   {
      v->S[n] = default_S[n] + user_S[n];
   }

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
   grid_ local_grid;

   local_grid.x[0] = grid.time;

#if DIM == 1

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = M_PI_2;

#elif DIM == 2

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = grid.X2[I[1]];

#elif DIM == 3 

   local_grid.x[1] = grid.X1[I[0]];
   local_grid.x[2] = grid.X2[I[1]];
   local_grid.x[3] = grid.X3[I[2]];
   
#endif

   if(pm == 1)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1p;
            local_grid.x[1] = grid.X1p[I[0]];
         break;

         case 'g':
            u = l->ux2p;
            local_grid.x[2] = grid.X2p[I[1]];
         break;

         case 'h':
            u = l->ux3p;
            local_grid.x[3] = grid.X3p[I[2]];
         break;
      }
   }
   if(pm == 0)
   {
      switch(flux)
      {
         case 'f':
            u = l->ux1m;
            local_grid.x[1] = grid.X1m[I[0]];
         break;

         case 'g':
            u = l->ux2m;
            local_grid.x[2] = grid.X2m[I[1]];
         break;

         case 'h':
            u = l->ux3m;
            local_grid.x[3] = grid.X3m[I[2]];
         break;
      }
   }

#if PHYSICS == RHD
   Get_Metric_Components(&local_grid);
#endif

   for(n = 0; n < eq; n++)
   {
      f->up[n] = u[1*eq + n];
      f->um[n] = u[0*eq + n];
   }

   Prim2Cons(f->qp,f->up,local_grid);
   Prim2Cons(f->qm,f->um,local_grid);

   switch(flux)
   {
      case 'f':
         Prim2FluxF(f->fp,dp,f->up,local_grid);

         Prim2FluxF(f->fm,dm,f->um,local_grid);
      break;

      case 'g':
         Prim2FluxG(f->fp,dp,f->up,local_grid);

         Prim2FluxG(f->fm,dm,f->um,local_grid);
      break;

      case 'h':
         Prim2FluxH(f->fp,dp,f->up,local_grid);

         Prim2FluxH(f->fm,dm,f->um,local_grid);
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
