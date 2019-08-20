/*
 *  aztekas initial module
 *  Date of creation: 11-07-2019 12:50:45
 *  author: Alejandro Aguayo Ortiz 
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   //////////////////////////////////
   // Kelvin-Helmholtz Instability //
   //////////////////////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(fabs(grid.X2[j]) >= x_0)
         {
            U(0,i,j) = nl;
            U(1,i,j) = pl;
            U(2,i,j) = vx1l*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j])); 
            U(3,i,j) = vx2l*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j]));
         }
         else if(fabs(grid.X2[j]) < x_0) 
         {
            U(0,i,j) = nr;
            U(1,i,j) = pr;
            U(2,i,j) = vx1r*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j])); 
            U(3,i,j) = vx2r*(1 + 0.01*cos(10*M_PI*grid.X1[i])*cos(10*M_PI*grid.X2[j]));
         }
      }
   }
   //////////////////////////////////
}
