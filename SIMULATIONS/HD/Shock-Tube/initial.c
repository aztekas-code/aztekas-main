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

void Initial()
{
   int n, i, j, k, cell;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      if(grid.X1[i] < x_0)
      {
         U(0,i) = nl;
         U(1,i) = pl;
         U(2,i) = vx1l;
      }
      else
      {
         U(0,i) = nr;
         U(1,i) = pr;
         U(2,i) = vx1r;
      }
   }

#elif DIM == 2

   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////

   #if INTERFACE == HORIZONTAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] < x_0)
         {
            U(0,i,j) = nl;
            U(1,i,j) = pl;
            U(2,i,j) = vx1l;
            U(3,i,j) = vx2l;
         }
         else 
         {
            U(0,i,j) = nr;
            U(1,i,j) = pr;
            U(2,i,j) = vx1r;
            U(3,i,j) = vx2r;
         }
      }
   }
   #elif INTERFACE == VERTICAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X2[j] > x_0)
         {
            U(0,i,j) = nl;
            U(1,i,j) = pl;
            U(2,i,j) = vx1l;
            U(3,i,j) = vx2l;
         }
         else 
         {
            U(0,i,j) = nr;
            U(1,i,j) = pr;
            U(2,i,j) = vx1r;
            U(3,i,j) = vx2r;
         }
      }
   }
   #elif INTERFACE == DIAGONAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X1[i] + grid.X2[j] - 1 < 0.0)
         {
            U(0,i,j) = nl;
            U(1,i,j) = pl;
            U(2,i,j) = vx1l;
            U(3,i,j) = vx2l;
         }
         else 
         {
            U(0,i,j) = nr;
            U(1,i,j) = pr;
            U(2,i,j) = vx1r;
            U(3,i,j) = vx2r;
         }
      }
   }
   #endif
#endif
}
