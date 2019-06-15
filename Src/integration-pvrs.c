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
#include"vector.h"

int INTEGRATION()
{
   int n, i, j, k;

   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
#if DIM == 1 
   
   RK1D(U,U,U1,U2,1);
   BOUNDARIES(U1);

   RK1D(U1,U,U1,U2,2);
   BOUNDARIES(U2);
   U = U2;
   
#elif DIM == 2 

   RK2D(U,U,U1,U2,1);
   BOUNDARIES(U1);

   RK2D(U1,U,U1,U2,2);
   BOUNDARIES(U2);
   U = U2;
   
#elif DIM == 3 
   
   funct_U2Q(Q,U);
   RK3D(U,Q,Q1,Q2,1);
   funct_Q2U(U1,Q1);
   BOUNDARIES(U1);
   RK3D(U1,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   BOUNDARIES(U);
   
#endif

   return 0;
}
