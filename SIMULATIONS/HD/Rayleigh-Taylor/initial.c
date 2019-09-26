/* 
 *  aztekas initial module
 *  Date of creation: 02-01-2019 12:50:45
 *  author: Alejandro Aguayo-Ortiz 
 */
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(grid.X2[j] <= 0.0)
         {
            U(RHO,i,j) =  rhod;
            U(PRE,i,j) =  pd - U(RHO,i,j)*grid.X2[j];
            U(VX1,i,j) =  vx1d;
            U(VX2,i,j) = -vx2d*(1 + cos(2*M_PI*grid.X1[i]/0.5))*(1 + cos(2*M_PI*grid.X2[j]/1.5))/4;
         }
         if(grid.X2[j] > 0.0)
         {
            U(RHO,i,j) =  rhou;
            U(PRE,i,j) =  pu - U(RHO,i,j)*grid.X2[j];
            U(VX1,i,j) =  vx1u;
            U(VX2,i,j) = -vx2u*(1 + cos(2*M_PI*grid.X1[i]/0.5))*(1 + cos(2*M_PI*grid.X2[j]/1.5))/4;
         }
      }
   }
}
