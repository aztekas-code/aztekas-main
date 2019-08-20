#include"main.h"
    
void Prim2Cons_All(double *q, double *u)
{
   int i, j, k;
   double rho, p, v_cov[3], v_con[3];
   double D, tau, S_cov[3], S_con[3];
   double Lorentz, W[3][3], U, VV, V[3];
   double P[eq+1];
   eos_ eos;
   gauge_ local_grid;

#if DIM == 1

   for(i = 0; i <= Nx1-0; i++)
   {
      local_grid.x[0] = grid.time;
      local_grid.x[1] = grid.X1[i];
      local_grid.x[2] = 0.0;
      local_grid.x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      local_grid.x[2] = M_PI_2;
      #endif
  
      Get_Metric_Components(&local_grid);

      rho      = u(0,i);
      p        = u(1,i);
      v_cov[0] = u(2,i);
      v_cov[1] = 0.0;
      v_cov[2] = 0.0;

      P[0] = rho;
      P[1] = p;

      v_con[0] = local_grid.gamma_con[0][0]*v_cov[0] + \
                 local_grid.gamma_con[0][1]*v_cov[1] + \
                 local_grid.gamma_con[0][2]*v_cov[2];
      v_con[1] = local_grid.gamma_con[1][0]*v_cov[0] + \
                 local_grid.gamma_con[1][1]*v_cov[1] + \
                 local_grid.gamma_con[1][2]*v_cov[2];
      v_con[2] = local_grid.gamma_con[2][0]*v_cov[0] + \
                 local_grid.gamma_con[2][1]*v_cov[1] + \
                 local_grid.gamma_con[2][2]*v_cov[2];
                 
      VV = v_cov[0]*v_con[0] + v_cov[1]*v_con[1] + v_cov[2]*v_con[2];

      Lorentz = 1.0/sqrt(1.0 - VV);

      EoS(&eos,P,local_grid);

      D        = rho*Lorentz;
      U        = rho*eos.h*Lorentz*Lorentz - p;
      tau      = U - D;
      S_cov[0] = rho*eos.h*Lorentz*Lorentz*v_cov[0];
      S_cov[1] = rho*eos.h*Lorentz*Lorentz*v_cov[1];
      S_cov[2] = rho*eos.h*Lorentz*Lorentz*v_cov[2];

      q(0,i) = D;
      q(1,i) = tau;
      q(2,i) = S_cov[0];
   }

#elif DIM == 2

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;
         #if POLAR == TRUE
         local_grid.x[2] = M_PI_2;
         #endif
     
         Get_Metric_Components(&local_grid);

         rho      = u(0,i,j);
         p        = u(1,i,j);
         v_cov[0] = u(2,i,j);
         v_cov[1] = u(3,i,j);
         v_cov[2] = 0.0;

         P[0] = rho;
         P[1] = p;

         v_con[0] = local_grid.gamma_con[0][0]*v_cov[0] + \
                    local_grid.gamma_con[0][1]*v_cov[1] + \
                    local_grid.gamma_con[0][2]*v_cov[2];
         v_con[1] = local_grid.gamma_con[1][0]*v_cov[0] + \
                    local_grid.gamma_con[1][1]*v_cov[1] + \
                    local_grid.gamma_con[1][2]*v_cov[2];
         v_con[2] = local_grid.gamma_con[2][0]*v_cov[0] + \
                    local_grid.gamma_con[2][1]*v_cov[1] + \
                    local_grid.gamma_con[2][2]*v_cov[2];
                    
         VV = v_cov[0]*v_con[0] + v_cov[1]*v_con[1] + v_cov[2]*v_con[2];

         Lorentz = 1.0/sqrt(1.0 - VV);

         EoS(&eos,P,local_grid);

         D        = rho*Lorentz;
         U        = rho*eos.h*Lorentz*Lorentz - p;
         tau      = U - D;
         S_cov[0] = rho*eos.h*Lorentz*Lorentz*v_cov[0];
         S_cov[1] = rho*eos.h*Lorentz*Lorentz*v_cov[1];
         S_cov[2] = rho*eos.h*Lorentz*Lorentz*v_cov[2];

         q(0,i,j) = D;
         q(1,i,j) = tau;
         q(2,i,j) = S_cov[0];
         q(3,i,j) = S_cov[1];
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1-0; i++)
   {
      for(j = 0; j <= Nx2-0; j++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;
       
         rho = u(0,i,j);
         p   = u(1,i,j);
         vx1 = u(2,i,j);
         vx2 = u(3,i,j);
         vx3 = u(4,i,j);

         P[0] = rho;
         P[1] = p;

         EoS(&eos,P,local_grid);

         E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
 
         q(0,i,j) = rho;
         q(1,i,j) = E;
         q(2,i,j) = rho*vx1;
         q(3,i,j) = rho*vx2;
         q(4,i,j) = rho*vx3;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            local_grid.x[0] = grid.time;
            local_grid.x[1] = grid.X1[i];
            local_grid.x[2] = grid.X2[j];
            local_grid.x[3] = grid.X3[k];
        
            rho = u(0,i,j,k);
            p   = u(1,i,j,k);
            vx1 = u(2,i,j,k);
            vx2 = u(3,i,j,k);
            vx3 = u(4,i,j,k);

            P[0] = rho;
            P[1] = p;

            EoS(&eos,P,local_grid);

            E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
          
            q(0,i,j,k) = rho;
            q(1,i,j,k) = E;
            q(2,i,j,k) = rho*vx1;
            q(3,i,j,k) = rho*vx2;
            q(4,i,j,k) = rho*vx3;
         }
      }
   }

#endif
}
