#include"main.h"

void User_Boundaries(double *B)
{
   int i, j, k, n, cell;
   gauge_ local_grid;

#if DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];

         Get_Metric_Components(&local_grid);
    #if POLAR == FALSE
         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
    #elif POLAR == TRUE
         if(grid.X2[j] >= M_PI_2 && grid.X2[j] <= 3.0*M_PI_2 && i >= Nx1-gc)
    #endif
         {
            B(0,i,j) =  density_inf;
            B(1,i,j) =  (K - 1.0)*B(0,i,j)*pow(velocity_inf/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_inf/Mach,2.0));
            B(2,i,j) =  velocity_inf*cos(grid.X2[j])/sqrt(local_grid.gamma_con[0][0]);
            B(3,i,j) = -velocity_inf*sin(grid.X2[j])/sqrt(local_grid.gamma_con[1][1]);
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
            B(0,i,j) =  density_inf;
            B(1,i,j) =  (K - 1.0)*B(0,i,j)*pow(velocity_inf/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_inf/Mach,2.0));
            B(2,i,j) =  velocity_inf*cos(grid.X2[j])/sqrt(local_grid.gamma_con[0][0]);
            B(3,i,j) = -velocity_inf*sin(grid.X2[j])/sqrt(local_grid.gamma_con[1][1]);
            B(4,i,j) = 0.0;
         }
      }
   }

#endif
}
