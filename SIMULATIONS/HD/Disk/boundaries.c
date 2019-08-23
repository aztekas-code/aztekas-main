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
   Periodic(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = r_dot_0;
            B(3,i,j) = phi_dot_0/grid.X1[i];
         }
      }
   }

   return 0;
}
