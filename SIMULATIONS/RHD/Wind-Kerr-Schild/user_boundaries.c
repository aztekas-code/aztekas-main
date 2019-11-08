#include"main.h"

void User_Boundaries(double *B)
{
   int i, j, k, n, cell;
   gauge_ local_grid;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         local_grid.x[1] = grid.X1[i];
         #if POLAR == FALSE
         local_grid.x[2] = grid.X2[j];
         #elif POLAR == TRUE
         local_grid.x[2] = M_PI_2;
         #endif

         double r = grid.X1[i];
         double p = grid.X2[j];
         double t = M_PI_2;
         double M = Black_Hole_Mass;
         double a = Black_Hole_Spin;

         Get_Metric_Components(&local_grid);

         double rho2 = r*r + a*a*cos(t)*cos(t);
         
         double gamrr = 1.0 + 2.0*M*r/rho2;
         double gamrp = -a*(gamrr)*sin(t)*sin(t);
         double gampp = sin(t)*sin(t)*(rho2 + a*a*(gamrr)*sin(t)*sin(t));

         double F1 = 1.0/sqrt(gamrr);
         double F4 = -2.0*gamrp/(sqrt(gamrr)*gampp);
         double F3 = (F1*gamrr + F4*gamrp)/sqrt((gamrr*gampp - -gamrp*gamrp));
         double F2 = (F4*gampp + F1*gamrp)/sqrt((gamrr*gampp - -gamrp*gamrp));

         double vr =  F1*velocity_0*cos(grid.X2[j]) + F2*velocity_0*sin(grid.X2[j]);
         double vp = -F3*velocity_0*sin(grid.X2[j]) + F4*velocity_0*cos(grid.X2[j]);

   #if POLAR == FALSE
         if(grid.X2[j] >= M_PI_2 && i >= Nx1-gc)
   #elif POLAR == TRUE
         if(grid.X2[j] >= M_PI_2 && grid.X2[j] <= 3.0*M_PI_2 && i >= Nx1-gc)
   #endif
         {
            B(0,i,j) =  density_0;
            B(1,i,j) =  (K - 1.0)*U(0,i,j)*pow(velocity_0/Mach,2.0)/(K*(K - 1.0) - K*pow(velocity_0/Mach,2.0));
            B(2,i,j) = gamrr*vr + gamrp*vp; 
            B(3,i,j) = gamrp*vr + gampp*vp;
         }
      }
   }
}
