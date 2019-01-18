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

void OUTFLOW(double *B)
{
   int i, j, k, n, cell;

   if(dim == 1)
   {
      for(n = 0; n < eq; n++)
      {
         for(cell = 0; cell < gc; cell++)
         {
            if(alfa > 0)
            {
               B[c1(n,gc)] = B[c1(n,gc+1)];
            }

            B[c1(n,cell)] = B[c1(n,gc)];
            B[c1(n,Nx1-cell)] = B[c1(n,Nx1-gc)];
         }
      }
   }
   else if(dim == 2)
   {
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1; i++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               if(alfa > 0)                                                       
               {                                                                   
                  B[c2(n,i,gc)] = B[c2(n,i,gc+1)];                                 
                  B[c2(n,i,Nx2-gc)] = B[c2(n,i,Nx2-gc-1)];                                 
               }                     

               B[c2(n,i,cell)] = B[c2(n,i,gc)];
               B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-gc)];
            }
         }

         for(j = 0; j <= Nx2; j++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               if(alfa > 0)
               {
                  B[c2(n,gc,j)] = B[c2(n,gc+1,j)];
                  B[c2(n,Nx1-gc,j)] = B[c2(n,Nx1-gc-1,j)];
               }

               B[c2(n,cell,j)] = B[c2(n,gc,j)];
               B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-gc,j)];
            }
         }
      }
   }
}


void PERIODIC(double *B, int r, int l, int u, int d, int f, int b)
{
   int i, j, k, n, cell;
   if(dim == 1)
   {
      for(n = 0; n < eq; n++)
      {
         for(i = 0; i <= Nx1-0; i++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               if(l == 1)
               {
                  B[c1(n,cell)] = B[c1(n,Nx1-2*gc+cell+1)];
               }
               if(r == 1)
               {
                  B[c1(n,Nx1-cell)] = B[c1(n,2*gc-cell-1)];
               }
            }
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
               if(d == 1)
               {
                  B[c2(n,i,cell)] = B[c2(n,i,Nx2-2*gc+cell+1)];
               }
               if(u == 1)
               {
                  B[c2(n,i,Nx2-cell)] = B[c2(n,i,2*gc-cell-1)];
               }
            }
         }

         for(j = 0; j <= Nx2; j++)
         {
            for(cell = 0; cell < gc; cell++)
            {
               if(l == 1)
               {
                  B[c2(n,cell,j)] = B[c2(n,Nx1-2*gc+cell+1,j)];
               }
               if(r == 1)
               {
                  B[c2(n,Nx1-cell,j)] = B[c2(n,2*gc-cell-1,j)];
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////
// REFLECTION: 
//    This function stablish a reflective condition in the x1, -x1, 
//    x2, -x2, x3, -x3 boundary, respectibly.
//
//    The conditions fix the point gc and Nxi-gc, i.e., the first 
//    and last points of the integration domain, to fulfil the 
//    condition: Tangent velocity = 0.0 and V[gc-1] = - V[gc+1]
////////////////////////////////////////////////////////////////////

void REFLECTION(double *B, int r, int l, int u, int d, int f, int b)
{
   int i, j, k, n, cell;
   if(dim == 1)
   {
      for(i = 0; i <= Nx1-0; i++)
      {
         for(cell = 0; cell < gc; cell++)
         {
            if(l == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c1(n,cell)] = B[c1(n,2*gc-cell-1)];
               }

               B[c1(2,cell)] = -B[c1(2,2*gc-cell-1)];
            }
            if(r == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c1(n,Nx1-cell)] = B[c1(n,Nx1-2*gc+cell+1)];
               }

               B[c1(2,Nx1-cell)] = -B[c1(2,Nx1-2*gc+cell+1)];
            }
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1-0; i++)
      {
         for(cell = 0; cell < gc; cell++)
         {
            if(d == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c2(n,i,cell)] = B[c2(n,i,2*gc-cell)];
               }

               B[c2(3,i,cell)] = -B[c2(3,i,2*gc-cell)];
               B[c2(3,i,gc)] = 0.0;
            }

            if(u == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-2*gc+cell)];
               }

               B[c2(3,i,Nx2-cell)] = -B[c2(3,i,Nx2-2*gc+cell)];
               B[c2(3,i,Nx2-gc)] = 0.0;
            }
         }
      }

      for(j = 0; j <= Nx2; j++)
      {
         for(cell = 0; cell < gc; cell++)
         {
            if(l == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c2(n,cell,j)] = B[c2(n,2*gc-cell-1,j)];
               }

               B[c2(2,cell,j)] = -B[c2(2,2*gc-cell-1,j)];
               B[c2(2,gc,j)] = 0.0;
            }

            if(r == 1)
            {
               for(n = 0; n < eq; n++)
               {
                  B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-2*gc+cell,j)];
               }

               B[c2(2,Nx1-cell,j)] = -B[c2(2,Nx1-2*gc+cell,j)];
               B[c2(2,Nx1-gc,j)] = 0.0;
            }
         }
      }
   }
}

void JET_LAUNCH(double *B)
{
   int i, j, k, n;
   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         if(fabs(X1[i]) <= r_jet)
         {
            B[c1(0,i)] = n_jet;
            B[c1(1,i)] = p_jet;
            B[c1(2,i)] = vx1_jet;
         }
      }
   }
   else if(dim == 2)
   {
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
   }
}

void IN_OUT_BOUND(double *B)
{
   int i, j, k, n, cell;
   double x1, x2, x3, r;

   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         r = X1[i];
         r_in = X1[gc+2];

         if(r > r_out)
         {
            B[c1(0,i)] = density_0;
            B[c1(1,i)] = pressure_0;
            B[c1(2,i)] = velocity_0;
         }
         else if(r <= r_in)
         {
            B[c1(0,i)] = B[c1(0,gc+2)];//density_0;
            B[c1(1,i)] = B[c1(1,gc+2)];//pressure_0;
            B[c1(2,i)] = B[c1(2,gc+2)];//0.0;
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {   
         for(j = 0; j <= Nx2; j++)
         {
            x1 = X1[i];
            x2 = X2[j];

            // Cart or Cyl
            if (alfa <= 1)
            {
               r = sqrt(x1*x1 + x2*x2);
               r_in = X1[gc+2];

               if(r >= r_out)
               {
                  B[c2(0,i,j)] = density_0;
                  B[c2(1,i,j)] = pressure_0;
                  B[c2(2,i,j)] = velocity_0*(x1/r);
                  B[c2(3,i,j)] = velocity_0*(x2/r);
               }
               else if(r <= r_in)
               {
                  B[c2(0,i,j)] = density_0;
                  B[c2(1,i,j)] = pressure_0;
                  B[c2(2,i,j)] = 0.0;
                  B[c2(3,i,j)] = 0.0;
               }
            }
            else if(alfa == 2)
            {
               r = x1;                                                          
               r_in = X1[gc+2]; 
                                                                           
               if(r >= r_out)                                                    
               {                                                                
                  B[c2(0,i,j)] = density_0;                                     
                  B[c2(1,i,j)] = pressure_0;                                    
                  B[c2(2,i,j)] = velocity_0;                                    
                  B[c2(3,i,j)] = 0.0;                                           
               }                                                                
               else if(r <= r_in)                                               
               {                                                                
                  B[c2(0,i,j)] = density_0;                                     
                  B[c2(1,i,j)] = pressure_0;                                    
                  B[c2(2,i,j)] =  0.0;                                          
                  B[c2(3,i,j)] =  0.0;    
               }
            }
         }
      }
   }
}

void WIND_BOUND(double *B)
{
   int i, j, k, cell;
   double r;

   if (dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            x1 = X1[i];
            x2 = X2[j];

            if(alfa <= 1)
            {
               r = sqrt(x1*x1 + x2*x2);

               if(X2[j] <= X2[2*gc])
               {
                  B[c2(0,i,j)] = density_0;                                     
                  B[c2(1,i,j)] = pressure_0;                                    
                  B[c2(2,i,j)] = 0.0;                                           
                  B[c2(3,i,j)] = velocity_0;                                    
               }
               else if(r <= r_in)
               {
                  B[c2(0,i,j)] = density_0;
                  B[c2(1,i,j)] = pressure_0;
                  B[c2(2,i,j)] = 0;
                  B[c2(3,i,j)] = 0;
               }
            }
            else if(alfa == 2)
            {
               r = x1;

               if(X2[j] >= M_PI_2 && r >= r_out)
               {
                  B[c2(0,i,j)] = density_0;                                     
                  B[c2(1,i,j)] = pressure_0;                                    
                  B[c2(2,i,j)] =  velocity_0*cos(X2[j])*(sqrt(1 + 2/X1[i]));                                           
                  B[c2(3,i,j)] = -velocity_0*sin(X2[j])*X1[i];                                    
               }
               else if(r <= r_in)
               {
                  B[c2(0,i,j)] = density_0;
                  B[c2(1,i,j)] = pressure_0;
                  B[c2(2,i,j)] = 0.0;
                  B[c2(3,i,j)] = 0.0;
               }
            }
         }
      }
   }
}
