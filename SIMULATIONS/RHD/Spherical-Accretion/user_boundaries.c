/* 
 *  aztekas user boundaries module
 *  Date of creation: 27-11-2019 10:44:40
 *  author: Alejandro Aguayo-Ortiz 
 */

#include"main.h"

void User_Boundaries(double *B)
{
   double r;
   double t;
   double R, z;
   double alpha, betar, gamrr;
   double vr;
   double test, eta;
   gauge_ local_grid;

#if DIM == 1

   for(int i = 0; i <= Nx1; i++)
   {
      r = grid.X1[i];

      local_grid.x[1] = r;
      local_grid.x[2] = M_PI_2;

      #if PHYSICS == RHD
      Get_Metric_Components(&local_grid);
      #endif 

      gamrr = local_grid.gamma_con[0][0];
      betar = local_grid.beta_con[0];
      alpha = local_grid.lapse;

      vr = gamrr*B(VX1,i) - betar/alpha;

      if(i >= Nx1-gc)
      {
         B(RHO,i) = 1.0;
         //B(PRE,j) = ((K - 1.0)*pow(B(RHO,i,j),K)*cs*cs)/(K*(K - 1.0) - K*cs*cs)/pow(1.0/density_0,K-1.0);
         B(PRE,i) = Temp*B(RHO,i)*pow(B(RHO,i)/(1.0/density_0),K-1.0);
         //B(VX1,i) = betar/alpha/gamrr;
         B(VX1,i) = velocity_0;
      }
   }

#elif DIM == 2 || DIM == 4
   for(int j = 0; j <= Nx2; j++)
    {
      for(int i = 0; i <= Nx1; i++)
       {
         r = grid.X1[i];
         t = grid.X2[j];

         local_grid.x[1] = r;
         local_grid.x[2] = t;

         #if PHYSICS == RHD
         Get_Metric_Components(&local_grid);
         #endif 

         gamrr = local_grid.gamma_con[0][0];
         betar = local_grid.beta_con[0];
         alpha = local_grid.lapse;

         vr = gamrr*B(VX1,i,j) - betar/alpha;

         if(i >= Nx1-gc)
          {
            B(RHO,i,j) = 1.0;
            //B(PRE,i,j) = ((K - 1.0)*pow(B(RHO,i,j),K)*cs*cs)/(K*(K - 1.0) - K*cs*cs)/pow(1.0/density_0,K-1.0);
            B(PRE,i,j) = Temp*B(RHO,i,j)*pow(density_0,K-1.0);
            //B(VX1,i,j) = betar/alpha/gamrr;
            B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = 0.0;
         }
      }
   }
#endif
}
