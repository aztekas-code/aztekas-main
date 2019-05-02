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
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

int MESH()
{
   int i, j, k;

#if dim == 1
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
#elif dim == 2 || dim == 4
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
#elif dim == 3
   dx1 = (x1max - x1min)/((double)Nx1-2*gc);
   dx2 = (x2max - x2min)/((double)Nx2-2*gc);
   dx3 = (x3max - x3min)/((double)Nx3-2*gc);
#endif

#if dim == 1
   
   for(i = 0; i <= Nx1; i++)
   {
      X1[i] = x1min + (i-gc)*(dx1);
      X1p[i] = x1min + (i+0.5-gc)*(dx1);
      X1m[i] = x1min + (i-0.5-gc)*(dx1);

      if(logmesh == 1)
      {
         X1[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2*gc)) - 1;
         X1p[i] = x1min + exp(log((x1max - x1min + 1.0))*(i+0.5-gc)/(Nx1-2*gc)) - 1;
         X1m[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-0.5-gc)/(Nx1-2*gc)) - 1;
      }
   }
   
#elif dim == 2  || dim == 4
   
   for(i = 0; i <= Nx1; i++)
   {
      X1[i]  = x1min + (i-gc)*(dx1);
      X1p[i] = x1min + (i+0.5-gc)*(dx1);
      X1m[i] = x1min + (i-0.5-gc)*(dx1);

      if(logmesh == 1)
      {
         X1[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2*gc)) - 1;
         X1p[i] = x1min + exp(log((x1max - x1min + 1.0))*(i+0.5-gc)/(Nx1-2*gc)) - 1;
         X1m[i] = x1min + exp(log((x1max - x1min + 1.0))*(i-0.5-gc)/(Nx1-2*gc)) - 1;
      }
   }

   for(j = 0; j <= Nx2; j++)
   {
      X2[j]  = x2min + (j-gc)*(dx2);
      X2p[j] = x2min + (j+0.5-gc)*(dx2);
      X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }
   
#elif dim == 3 
   
   for(i = 0; i <= Nx1; i++)
   {
      X1[i]  = x1min + (i-gc)*(dx1);
      X1p[i] = x1min + (i+0.5-gc)*(dx1);
      X1m[i] = x1min + (i-0.5-gc)*(dx1);
   }

   for(j = 0; j <= Nx2; j++)
   {
      X2[j]  = x2min + (j-gc)*(dx2);
      X2p[j] = x2min + (j+0.5-gc)*(dx2);
      X2m[j] = x2min + (j-0.5-gc)*(dx2);
   }

   for(k = 0; k <= Nx3; k++)
   {
      X3[k]  = x3min + (k-gc)*(dx3);
      X3p[k] = x3min + (k+0.5-gc)*(dx3);
      X3m[k] = x3min + (k-0.5-gc)*(dx3);
   }

#endif
   return 0;
}
