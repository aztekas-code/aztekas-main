#include"main.h"
    
void Prim2FluxF(double *f, double *v, double *u, gauge_ local_grid)
{
   int i;
   double rho, p, v_cov[3], v_con[3];
   double D, U, tau, Lorentz, h, cs, VV;
   double W[3][3], S_cov[3], S_con[3], V[3];
   double gamma, beta, alpha, vel;
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

   v_con[0] = local_grid.gamma_con[0][0]*v_cov[0] + \
              local_grid.gamma_con[0][1]*v_cov[1] + \
              local_grid.gamma_con[0][2]*v_cov[2];
   v_con[1] = local_grid.gamma_con[1][0]*v_cov[0] + \
              local_grid.gamma_con[1][1]*v_cov[1] + \
              local_grid.gamma_con[1][2]*v_cov[2];
   v_con[2] = local_grid.gamma_con[2][0]*v_cov[0] + \
              local_grid.gamma_con[2][1]*v_cov[1] + \
              local_grid.gamma_con[2][2]*v_cov[2];

   VV = v_con[0]*v_cov[0] + v_con[1]*v_cov[1] + v_con[2]*v_cov[2];

   Lorentz = 1/sqrt(1 - VV);

   EoS(&eos,u,local_grid);

   h  = eos.h;
   cs = eos.cs;

   D   = rho*Lorentz;
   U   = rho*h*Lorentz*Lorentz - p;
   tau = U - D;

   for(i = 0; i < 3; i++)
   {
      V[i] = local_grid.lapse*v_con[i] - local_grid.beta_con[i];
      S_cov[i] = rho*eos.h*Lorentz*Lorentz*v_cov[i];
      S_con[i] = rho*eos.h*Lorentz*Lorentz*v_con[i];
   }

   W[0][0] = S_con[0]*v_cov[0] + p;
   W[0][1] = S_con[0]*v_cov[1];
   W[0][2] = S_con[0]*v_cov[2];

   f[0] = D*V[0];
   f[1] = local_grid.lapse*(S_con[0] - v_con[0]*D) - local_grid.beta_con[0]*tau;
   f[2] = W[0][0];
   f[3] = pow(Lorentz,2.0)*h*rho*v_cov[0]*v_cov[1];
   f[4] = pow(Lorentz,2.0)*h*rho*v_cov[0]*v_cov[2];

   gamma = local_grid.gamma_con[0][0];
   beta  = local_grid.beta_con[0];
   alpha = local_grid.lapse;
   vel   = v_con[0];

   v[0] = -(fabs(cs)*sqrt(((pow(VV,2.0)-VV)*pow(cs,2.0)-VV+1)*gamma+(1-VV)*pow(vel,2.0)*pow(cs,2.0    )+(VV-1)*pow(vel,2.0))+(beta*VV-vel)*pow(cs,2.0)+vel-beta)/(VV*pow(cs,2.0)-1);
   v[1] = (fabs(cs)*sqrt(((pow(VV,2.0)-VV)*pow(cs,2.0)-VV+1)*gamma+(1-VV)*pow(vel,2.0)*pow(cs,2.0)    +(VV-1)*pow(vel,2.0))+(vel-beta*VV)*pow(cs,2.0)-vel+beta)/(VV*pow(cs,2.0)-1);
   v[2] = alpha*vel-beta;
}
