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

#if dim == 1
   OUTFLOW(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(n = 0; n < eq; n++)
      {
         roundgen(&B[c1(n,i)]);
      }
   }
#elif dim == 2
   OUTFLOW(B);
   REFLECTIVE(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(n = 0; n < eq; n++)
         {
            roundgen(&B[c2(n,i,j)]);
         }
      }
   }
#endif

   return 0;
}
