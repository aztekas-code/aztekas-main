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
#include<string.h>
#include"main.h"
#include"param.h"

void INITIAL(double *dtprint)
{
   int n, i, j, k, cell;

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

#if dim == 1
   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      if(fabs(X1[i]) < x_0)
      {
         U[c1(0,i)] = nl;
         U[c1(1,i)] = pl;
         U[c1(2,i)] = vx1l;
      }
      else 
      {
         U[c1(0,i)] = nr;
         U[c1(1,i)] = pr;
         U[c1(2,i)] = vx1r;
      }
   }

#elif dim == 2
   /////////////////////////////
   //-------Riemann-2D--------//
   /////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(X1[i] < x_0)
         {
            U[c2(0,i,j)] = nl;
            U[c2(1,i,j)] = pl;
            U[c2(2,i,j)] = vx1l;
            U[c2(3,i,j)] = vx2l;
         }
         else 
         {
            U[c2(0,i,j)] = nr;
            U[c2(1,i,j)] = pr;
            U[c2(2,i,j)] = vx1r;
            U[c2(3,i,j)] = vx2r;
         }
      }
   }

#endif
}
