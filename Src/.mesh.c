/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version gc of the License, or
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

int Mesh()
{
   int i, j, k;

#if DIM == 1
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
#elif DIM == 2 || DIM == 4
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
#elif DIM == 3
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
   dx3 = (x3max - x3min)/((double)Nx3-2*gc);
#endif

#if DIM == 1
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i] = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);

      if(logmesh == 1)
      {
         grid.X1[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2*gc)) - 1;
         grid.X1p[i] = x1min + exp(log((x1max - x1min + 1.0))*(i+0.5-gc)/(Nx1-2*gc)) - 1;
         grid.X1m[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-0.5-gc)/(Nx1-2*gc)) - 1;
      }
   }
   
#elif DIM == 2  || DIM == 4
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i]  = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);

      if(logmesh == 1)
      {
         grid.X1[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2*gc)) - 1;
         grid.X1p[i] = x1min + exp(log((x1max - x1min + 1.0))*(i+0.5-gc)/(Nx1-2*gc)) - 1;
         grid.X1m[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-0.5-gc)/(Nx1-2*gc)) - 1;
      }
   }

   for(j = 0; j <= Nx2; j++)
   {
      grid.X2[j]  = x2min + (j-gc)*(dx2);
      grid.X2p[j] = x2min + (j+0.5-gc)*(dx2);
      grid.X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }
   
#elif DIM == 3 
   
   for(i = 0; i <= Nx1; i++)
   {
      grid.X1[i]  = x1min + (i-gc)*(dx1);
      grid.X1p[i] = x1min + (i+0.5-gc)*(dx1);
      grid.X1m[i] = x1min + (i-0.5-gc)*(dx1);
   }

   for(j = 0; j <= Nx2; j++)
   {
      grid.X2[j]  = x2min + (j-gc)*(dx2);
      grid.X2p[j] = x2min + (j+0.5-gc)*(dx2);
      grid.X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }

   for(k = 0; k <= Nx3; k++)
   {
      grid.X3[k]  = x3min + (k-gc)*(dx3);
      grid.X3p[k] = x3min + (k+0.5-gc)*(dx3);
      grid.X3m[k] = x3min + (k-0.5-gc)*(dx3);
   }

#endif
   return 0;
}
