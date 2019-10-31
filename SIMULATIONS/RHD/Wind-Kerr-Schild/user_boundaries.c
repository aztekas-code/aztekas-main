#include"main.h"

void User_Boundaries(double *B)
{
   gauge_ local_grid;

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];

         Get_Metric_Components(&local_grid);

         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
         {
            B(RHO,i,j) =  density_0;
            B(PRE,i,j) =  (K - 1.0)*B(RHO,i,j)*pow(velocity_0/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_0/Mach,2.0));
            B(VX1,i,j) =  velocity_0*cos(grid.X2[j])/sqrt(local_grid.gamma_con[0][0]);
            B(VX2,i,j) = -velocity_0*sin(grid.X2[j])/sqrt(local_grid.gamma_con[1][1]);
         }
      }
   }
}
