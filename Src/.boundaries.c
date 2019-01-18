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

int BOUNDARIES(double *B)
{
   int n, i, j, k, cell;

   // Set by default outflow boundaries
   OUTFLOW(B);

   switch(TEST)
   {
      //x1 (right), -x1 (left), x2 (up), -x2 (down), x3, -x3
      case 0: //Cutom boundaries
         BOUND_CUSTOM(B);
      break;

      case 1:
//         REFLECTION(B,0,0,1,1,0,0);
      break;

      case 2: //Kelvin-Helmholtz 
         PERIODIC(B,1,1,1,1,0,0);
      break;

      case 3: //Jet
         REFLECTION(B,0,1,0,0,0,0);
         JET_LAUNCH(B);
      break;

      case 4: //Spherical accretion
         // 0,1,0,1,0,0 for cyl 2D
         // 0,1,0,0,0,0 for sph 2D
         if(alfa <= 1)
         {
            REFLECTION(B,0,1,0,1,0,0);
         }
         else if(alfa == 2)
         {
            REFLECTION(B,0,0,1,1,0,0);
         }
         IN_OUT_BOUND(B);
      break;

      case 5:
         if(alfa <= 1)
         {
            REFLECTION(B,0,1,0,0,0,0);
         }
         else if(alfa == 2)
         {
            REFLECTION(B,0,0,1,1,0,0);
         }
         WIND_BOUND(B);
      break;
   }

/*
   if(dim == 1)
   {
      for(n = 0; n < eq; n++)
      {
         for(cell = 0; cell < gc; cell++)
         {
            B[c1(n,cell)] = B[c1(n,gc)];
            B[c1(n,Nx1-cell)] = B[c1(n,Nx1-gc)];
         }
      }
   }
   else if(dim == 2)
   {
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               B[c2(n,i,cell)] = B[c2(n,i,gc)];
               B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-gc)];
            }
         }

         for(j = 0; j <= Nx2; j++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               B[c2(n,cell,j)] = B[c2(n,gc,j)];
               B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-gc,j)];
            }
         }
      }

      ////////////////////////////////////////////////
      //-------------Kelvin-Helmholtz---------------//
      ////////////////////////////////////////////////

      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
            //Reflection
            for(cell = 0; cell < gc; cell++)
            {
               B[c2(n,i,cell)] = B[c2(n,i,Nx2-2*gc-cell-1)];
               B[c2(n,i,Nx2-cell)] = B[c2(n,i,2*gc-cell-1)];
            }
         }

         for(j = 0; j <= Nx2; j++)
         {
            B[c2(n,2,j)] = B[c2(n,Nx1-3,j)];
            B[c2(n,1,j)] = B[c2(n,Nx1-4,j)];
            B[c2(n,0,j)] = B[c2(n,Nx1-5,j)];

            B[c2(n,Nx1-2,j)] = B[c2(n,3,j)];
            B[c2(n,Nx1-1,j)] = B[c2(n,4,j)];
            B[c2(n,Nx1  ,j)] = B[c2(n,5,j)];
         }
      }

      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      ///*-------------ACRETION-------------------*///
      ////////////////////////////////////////////////
/*
      for(j = 0; j <= Nx2; j++)
      {
         B[c2(2,2,j)] = -B[c2(2,4,j)];
         B[c2(2,1,j)] = -B[c2(2,5,j)];
         B[c2(2,0,j)] = -B[c2(2,6,j)];

         B[c2(2,3,j)] = 0.5*(B[c2(2,2,j)] + B[c2(2,4,j)]);
      }

      //Viento
      for(i = 0; i <= Nx1; i++)
      {
         B[c2(0,i,5)] = 1.0;
         B[c2(1,i,5)] = cs*cs*B[c2(0,i,5)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,5)] = 0.0;
         B[c2(3,i,5)] = v;

         B[c2(0,i,4)] = 1.0;
         B[c2(1,i,4)] = cs*cs*B[c2(0,i,4)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,4)] = 0.0;
         B[c2(3,i,4)] = v;

         B[c2(0,i,3)] = 1.0;
         B[c2(1,i,3)] = cs*cs*B[c2(0,i,3)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,3)] = 0.0;
         B[c2(3,i,3)] = v;

         B[c2(0,i,2)] = 1.0;
         B[c2(1,i,2)] = cs*cs*B[c2(0,i,2)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,2)] = 0.0;
         B[c2(3,i,2)] = v;

         B[c2(0,i,1)] = 1.0;
         B[c2(1,i,1)] = cs*cs*B[c2(0,i,1)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,1)] = 0.0;
         B[c2(3,i,1)] = v;

         B[c2(0,i,0)] = 1.0;
         B[c2(1,i,0)] = cs*cs*B[c2(0,i,0)]*(K-1)/(K*(K-1) - cs*cs*K);
         B[c2(2,i,0)] = 0.0;
         B[c2(3,i,0)] = v;
      }

      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            double RR = sqrt(X1[i]*X1[i] + X2[j]*X2[j]);
            double R  = X1[i];
            double z  = X2[j];
            double r  = MM/(cs*cs + v*v);

            B[c2(0,i,j)] = fabs(B[c2(0,i,j)]);
            B[c2(1,i,j)] = fabs(B[c2(1,i,j)]);

            if(RR <= rr)
            {
               B[c2(0,i,j)] = 0.1;
               B[c2(1,i,j)] = cs*cs*B[c2(0,i,j)]*(K-1)/(K*(K-1) - cs*cs*K);
               B[c2(2,i,j)] = 0.0;
               B[c2(3,i,j)] = 0.0;
            }
         }
      }
*/
      ////////////////////////////////////////////////
      ///*---------JET-BOUNDARY-------------------*///
      ////////////////////////////////////////////////
/*
      for(j = 0; j <= Nx2; j++)
      {
         B[c2(2,3,j)] = -B[c2(2,4,j)];
         B[c2(2,2,j)] = -B[c2(2,5,j)];
         B[c2(2,1,j)] = -B[c2(2,6,j)];
         B[c2(2,0,j)] = -B[c2(2,7,j)];
      }

      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(fabs(X1[i]) <= 1.0)
            {
               if(fabs(X2[j]) <= 1.0)
               {
                  B[c2(0,i,j)] = n_jet;
                  B[c2(1,i,j)] = p_jet;
                  B[c2(2,i,j)] = vx1_jet;
                  B[c2(3,i,j)] = vx2_jet;
               }
            }
         }
      }
*/ 
/*   }
   else if(dim == 3)
   {
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
            for(k = 0; k <= Nx3-0; k++)
            {
               B[c3(n,i,2,k)] =  B[c3(n,i,3,k)];
               B[c3(n,i,1,k)] =  B[c3(n,i,4,k)];
               B[c3(n,i,0,k)] =  B[c3(n,i,5,k)];

               B[c3(n,i,Nx2-2,k)] = B[c3(n,i,Nx2-3,k)];
               B[c3(n,i,Nx2-1,k)] = B[c3(n,i,Nx2-4,k)];
               B[c3(n,i,Nx2  ,k)] = B[c3(n,i,Nx2-5,k)];
            }
         }

         for(j = 0; j <= Nx2-0; j++)
         {
            for(k = 0; k <= Nx3-0; k++)
            {
               B[c3(n,2,j,k)] = B[c3(n,3,j,k)];
               B[c3(n,1,j,k)] = B[c3(n,4,j,k)];
               B[c3(n,0,j,k)] = B[c3(n,5,j,k)];

               B[c3(n,Nx1-2,j,k)] = B[c3(n,Nx1-3,j,k)];
               B[c3(n,Nx1-1,j,k)] = B[c3(n,Nx1-4,j,k)];
               B[c3(n,Nx1  ,j,k)] = B[c3(n,Nx1-5,j,k)];
            }
         }

         for(i = 0; i <= Nx1-0; i++)
         {
            for(j = 0; j <= Nx2-0; j++)
            {
               B[c3(n,i,j,2)] = B[c3(n,i,j,3)];
               B[c3(n,i,j,1)] = B[c3(n,i,j,4)];
               B[c3(n,i,j,0)] = B[c3(n,i,j,5)];

               B[c3(n,i,j,Nx3-2)] = B[c3(n,i,j,Nx3-3)];
               B[c3(n,i,j,Nx3-1)] = B[c3(n,i,j,Nx3-4)];
               B[c3(n,i,j,Nx3  )] = B[c3(n,i,j,Nx3-5)];
            }
         }
      }
   }
*/
   return 0;
}
