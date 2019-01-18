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

int c1(int n, int i)
{
   int c;

   c = n*(Nx1+1) + i;
   return c;
}

int c2(int n, int i, int j)
{
   int c;

   c = n*(Nx1+1)*(Nx2+1) + i*(Nx2+1) + j;
   return c;
}

int c3(int n, int i, int j, int k)
{
   int c;

   c = n*(Nx1+1)*(Nx2+1)*(Nx3+1) + i*(Nx2+1)*(Nx3+1) + j*(Nx3+1) + k;
   return c;
}
