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
#include"param.h"

int BOUNDARIES(double *B)
{
   int i, j, k, n, cell;

   ////////////////
   // OUTFLOW 1D //
   ////////////////
#if dim == 1
   for(cell = 0; cell < gc; cell++)
   {
      for(n = 0; n < eq; n++)
      {
         // r = 0 boundary in cylindrical and spherical coordinates
         #if alfa >= 1
         if(x1min == 0.0)
            B[c1(n,gc)] = B[c1(n,gc+1)];
         #endif

         B[c1(n,Nx1-cell)] = B[c1(n,Nx1-gc)]; //x1max
         B[c1(n,cell)] = B[c1(n,gc)]; //x1min
      }
   }
   ////////////////

#elif dim == 2
   ////////////////
   // OUTFLOW-2D //
   ////////////////

   // Outflow in X1 //
   for(j = 0; j <= Nx2; j++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         // r = 0 boundary in cylindrical and spherical coordinates
         #if alfa >= 1
         if(x1min == 0.0) 
             B[c2(n,gc,j)] = B[c2(n,gc+1,j)];
         #endif

         for(n = 0; n < eq; n++)
         {
            B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-gc,j)]; //x1max
            B[c2(n,cell,j)] = B[c2(n,gc,j)]; //x1min
         }
      }
   }
   ///////////////////

   // Outflow in X2 //
   for(i = 0; i <= Nx1; i++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         // theta = 0 or M_PI boundary in spherical coordinates
         #if alfa == 2
         if(x2max/M_PI == 1) 
            B[c2(n,i,Nx2-gc)] = B[c2(n,i,Nx2-gc-1)];
         if(x2min/M_PI == 0) 
            B[c2(n,i,gc)] = B[c2(n,i,gc+1)];
         #endif

         for(n = 0; n < eq; n++)
         {
            B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-gc)]; //x2max
            B[c2(n,i,cell)] = B[c2(n,i,gc)]; //x2min

         }
      }
   }
#endif

   return 0;
}
