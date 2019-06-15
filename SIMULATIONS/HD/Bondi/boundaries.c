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

   Outflow(B);
   Reflection(B);

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B[c1(0,i)] = density_0;
         B[c1(1,i)] = pressure_0;
         B[c1(2,i)] = velocity_0;
      }
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B[c2(0,i,j)] = density_0;
            B[c2(1,i,j)] = pressure_0;
            B[c2(2,i,j)] = velocity_0;
            B[c2(3,i,j)] = 0.0;
         }
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B[c2(0,i,j)] = density_0;
            B[c2(1,i,j)] = pressure_0;
            B[c2(2,i,j)] = velocity_0;
            B[c2(3,i,j)] = 0.0;
            B[c2(4,i,j)] = 0.0;
         }
      }
   }

#endif

   return 0;
}
