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

int Boundaries(double *B)
{
   int i, j, k, n, cell;

#if DIM == 1
   Outflow(B);

#elif DIM == 2
   Outflow(B);
   
   // Linear extrapolation. Important for the diagonal shock tube
 
   #if INTERFACE == DIAGONAL
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(n = 0; n < eq; n++)
         {
            if(i > Nx1-gc)
            {
               B(n,i,j) = B(n,i-2,j) + ((grid.X1[i] - grid.X1[i-2])/(grid.X1[i-1] - grid.X1[i-2]))*(B(n,i-1,j) - B(n,i-2,j)); 
            }
            if(j > Nx2-gc)
            {
               B(n,i,j) = B(n,i,j-2) + ((grid.X2[j] - grid.X2[j-2])/(grid.X2[j-1] - grid.X2[j-2]))*(B(n,i,j-1) - B(n,i,j-2)); 
            }
         }
      }
   }

   for(i = Nx1; i >= 0; i--)
   {
      for(j = Nx2; j <= 0; j--)
      {
         for(n = 0; n < eq; n++)
         {
            if(i > Nx1-gc)
            {
               B(n,i,j) = B(n,i+2,j) + ((grid.X1[i] - grid.X1[i+2])/(grid.X1[i+1] - grid.X1[i+2]))*(B(n,i+1,j) - B(n,i+2,j)); 
            }
            if(j > Nx2-gc)
            {
               B(n,i,j) = B(n,i,j+2) + ((grid.X2[j] - grid.X2[j+2])/(grid.X2[j+1] - grid.X2[j+2]))*(B(n,i,j+1) - B(n,i,j+2)); 
            }
         }
      }
   }
   #endif 
   
#endif

   return 0;
}
