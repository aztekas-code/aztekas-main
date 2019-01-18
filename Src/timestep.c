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
#include"./Headers/main.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

double TIMESTEP()
{
   int i, j, k;
   double dtmin;
   double c, dt, cmax;
   double r;

   dtmin = 100000;

   if(dim == 1)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         c = sqrt(K*U[c1(1,i)] / (U[c1(0,i)]));
         dtmin = min(dx1/(fabs(U[c1(2,i)]) + fabs(c)),dtmin);
      }
   }
   else if(dim == 2)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            c = sqrt(K*U[c2(1,i,j)] / (U[c2(0,i,j)]));
            dtmin = min(dx1/(fabs(U[c2(2,i,j)]) + fabs(c)),dtmin);
            dtmin = min(dx2/(fabs(U[c2(3,i,j)]) + fabs(c)),dtmin);
         }
      }
   }
   else if(dim == 4)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            c = sqrt(K*U[c2(1,i,j)] / (U[c2(0,i,j)]));
            dtmin = min(dx1/(fabs(U[c2(2,i,j)]) + fabs(c)),dtmin);
            dtmin = min(dx1/(fabs(U[c2(4,i,j)]) + fabs(c)),dtmin);
            dtmin = min(dx2/(fabs(U[c2(3,i,j)]) + fabs(c)),dtmin);
            dtmin = min(dx2/(fabs(U[c2(4,i,j)]) + fabs(c)),dtmin);
         }
      }
   }
   else if(dim == 3)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            for(k = gc; k <= Nx3-gc; k++)
            {
               c = sqrt(K*U[c3(1,i,j,k)] / (U[c3(0,i,j,k)]));
               dtmin = min(dx1/(fabs(U[c3(2,i,j,k)]) + fabs(c)),dtmin);
               dtmin = min(dx2/(fabs(U[c3(3,i,j,k)]) + fabs(c)),dtmin);
               dtmin = min(dx3/(fabs(U[c3(4,i,j,k)]) + fabs(c)),dtmin);
            }
         }
      }
   }

#if PHYSICS == 1 //HD
   dt = cou*dtmin;
#elif PHYSICS == 2 //RHD
   #if dim == 1
   dt = cou*min(dx1,1000);
   #elif dim == 2 || dim == 4
   dt = cou*min(dx1,dx2);
   #endif 
#endif
   printf("%e\n",dt);

   return dt;
}
