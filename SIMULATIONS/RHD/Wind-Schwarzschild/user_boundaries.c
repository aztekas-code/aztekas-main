/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Relativistic Wind problem 
 * Creation Date : 05-12-2019
 * Last Modified : 17-12-2019 13:26:49
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
   int i, j, k, n, cell;
   gauge_ local_grid;

   Outflow(B);
#if POLAR == FALSE
   Reflection(B);
#elif POLAR == TRUE
   Periodic(B);
#endif

#if DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];

         Get_Metric_Components(&local_grid);

         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
         {
            B(0,i,j) =  density_0;
            B(1,i,j) =  (K - 1.0)*B(0,i,j)*pow(velocity_0/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_0/Mach,2.0));
            B(2,i,j) =  velocity_0*cos(grid.X2[j])/sqrt(local_grid.gamma_con[0][0]);
            B(3,i,j) = -velocity_0*sin(grid.X2[j])/sqrt(local_grid.gamma_con[1][1]);
         }
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];

         Get_Metric_Components(&local_grid);

         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
         {
            B(0,i,j) =  density_0;
            B(1,i,j) =  (K - 1.0)*B(0,i,j)*pow(velocity_0/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_0/Mach,2.0));
            B(2,i,j) =  velocity_0*cos(grid.X2[j])/sqrt(local_grid.gamma_con[0][0]);
            B(3,i,j) = -velocity_0*sin(grid.X2[j])/sqrt(local_grid.gamma_con[1][1]);
            B(4,i,j) = 0.0;
         }
      }
   }

#endif
}
