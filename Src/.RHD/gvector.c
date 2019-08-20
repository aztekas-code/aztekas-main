#include"main.h"
    
void Prim2FluxG(double *f, double *v, double *u, gauge_ local_grid)
{
   int i, j;
   double rho, p, v_cov[3], v_con[3];
   double D, tau, S_cov[3], S_con[3];
   double Lorentz, W[3][3], U, VV, V[3];
   eos_ eos;

   rho = u[0];
   p   = u[1];

#if DIM == 1
   v_cov[0] = u[2];
   v_cov[1] = 0.0;
   v_cov[2] = 0.0;
#elif DIM == 2
   v_cov[0] = u[2];
   v_cov[1] = u[3];
   v_cov[2] = 0.0;
#elif DIM == 3 || DIM == 4
   v_cov[0] = u[2];
   v_cov[1] = u[3];
   v_cov[2] = u[4];
#endif

   for(i = 0; i < 3; i++)
   {
      for(j = 0; j < 3; j++)
      {
         v_con[i] = local_grid.gamma_con[i][j]*v_cov[j];
         V[i] = local_grid.lapse*v_con[i] - local_grid.beta_con[i];
      }

      S_cov[i] = rho*eos.h*Lorentz*Lorentz*v_cov[i];
      S_con[i] = rho*eos.h*Lorentz*Lorentz*v_con[i];
   }

   Scalar_Contraction_Range1(&VV,v_cov,v_con);
   EoS(&eos,u,local_grid);

   Lorentz = 1.0/sqrt(1.0 - VV);

   D   = rho*Lorentz;
   U   = rho*eos.h*Lorentz*Lorentz - p;
   tau = U - D;

   W[1][0] = S_con[1]*v_con[0] + p*local_grid.gamma_con[1][0];
   W[1][1] = S_con[1]*v_con[1] + p*local_grid.gamma_con[1][1];
   W[1][2] = S_con[1]*v_con[2] + p*local_grid.gamma_con[1][2];

   f[0] = D*V[0];
   f[1] = local_grid.lapse*(S_con[1] - v_con[1]*D) - local_grid.beta_con[1]*tau;
   f[2] = local_grid.lapse*W[1][0] - local_grid.beta_con[1]*S_cov[0];
   f[3] = local_grid.lapse*W[1][1] - local_grid.beta_con[1]*S_cov[1];
   f[4] = local_grid.lapse*W[1][2] - local_grid.beta_con[1]*S_cov[2];

   double a2, vel, v2, gamma;
   double lambda_plus, lambda_minus;

   a2    = eos.cs*eos.cs;
   vel   = v_con[1];
   v2    = vel*vel;
   gamma = local_grid.gamma_con[1][1];

   lambda_plus = ((1 - a2)*vel + sqrt(a2*(1 - VV)*(1 - VV*a2)*gamma - (1 - a2)*v2))/(1 - VV*a2);
   lambda_plus = ((1 - a2)*vel - sqrt(a2*(1 - VV)*(1 - VV*a2)*gamma - (1 - a2)*v2))/(1 - VV*a2);

   v[0] = local_grid.lapse*vel - local_grid.beta_con[1];
   v[1] = local_grid.lapse*lambda_plus  - local_grid.beta_con[1];
   v[2] = local_grid.lapse*lambda_minus - local_grid.beta_con[1];
}
