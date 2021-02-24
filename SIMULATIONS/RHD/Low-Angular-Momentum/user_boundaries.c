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
   double alpha, beta, gamma;
   double vr;
   double test, eta;
   double M = Black_Hole_Mass;
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

      gamma = local_grid.gamma_con[0][0];
      beta = local_grid.beta_con[0];
      alpha = local_grid.lapse;

      vr = gamma*B(VX1,i) - beta/alpha;

      if(i >= Nx1-gc)
      {
         B(RHO,i) = 1.0;
         B(PRE,i) = Temp*B(RHO,i)*pow(B(RHO,i)/(1.0/density_0),K-1.0);
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

         gamma = local_grid.gamma_con[0][0];
         beta  = 2.0*M/r;//local_grid.beta_con[0];
         alpha = sqrt(r/(r + 2.0*M));//local_grid.lapse;
         vr = (r/(r + 2.0*M))*B(VX1,i,j);

         if(i >= Nx1-gc)
          {
            B(RHO,i,j) = 1.0;
            B(PRE,i,j) = Temp*B(RHO,i,j)*pow(density_0,K-1.0);
            //B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
            vr = (r/(r + 2.0*M))*B(VX1,i,j);
            B(VX3,i,j) = l_0*ftheta(t)*(alpha - beta*vr);
         }
      }
   }
#endif
}
