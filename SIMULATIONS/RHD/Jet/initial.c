/* 
 *  aztekas initial module
 *  Date of creation: 02-01-2019 18:32:00
 *  author: Alejandro Aguayo Ortiz 
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
         if(grid.X1[i] <= r_jet && grid.X2[j] <= z_jet)
         {
            U(0,i,j) = n_jet;
            U(1,i,j) = p_jet;
            U(2,i,j) = vx1_jet;
            U(3,i,j) = vx2_jet;
         }
         else
         {
            U(0,i,j) = n_atm;
            U(1,i,j) = p_atm;
            U(2,i,j) = vx1_atm;
            U(3,i,j) = vx2_atm;
         }
      }
   }
}
