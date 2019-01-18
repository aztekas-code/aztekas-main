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

int MxV(double *M, double *V, double *L)
{
   int n, m;
   double res=0.0;

   for(m = 0; m < eq; m++)
   {
      for(n = 0; n < eq; n++)
      {
         res += M[m*(eq) + n]*V[n];
      }

      L[m] = res;
      res = 0.0;
   }

   return 0;
}

void roundgen(double *num)
{
   double r;
   double bla;
   double decnum;
   int expnum;
   
   if(*num != 0.0e+00)
   {
      expnum = (int)floor(log(fabs(*num))/log(10.0));
      decnum = *num/pow(10,(double)expnum);
      decnum = roundf(decnum*1.0e+15)/1.0e+15;
      *num = decnum*pow(10,(double)expnum);
   }
}
