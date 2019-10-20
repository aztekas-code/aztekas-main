#include"main.h"
    
void Prim2Cons(double *q, double *u, gauge_ *local_grid)
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

   v_con[0] = local_grid->gamma_con[0][0]*v_cov[0] + \
              local_grid->gamma_con[0][1]*v_cov[1] + \
              local_grid->gamma_con[0][2]*v_cov[2];
   v_con[1] = local_grid->gamma_con[1][0]*v_cov[0] + \
              local_grid->gamma_con[1][1]*v_cov[1] + \
              local_grid->gamma_con[1][2]*v_cov[2];
   v_con[2] = local_grid->gamma_con[2][0]*v_cov[0] + \
              local_grid->gamma_con[2][1]*v_cov[1] + \
              local_grid->gamma_con[2][2]*v_cov[2];
              
   VV = v_cov[0]*v_con[0] + v_cov[1]*v_con[1] + v_cov[2]*v_con[2];

   Lorentz = 1.0/sqrt(1.0 - VV);

   EoS(&eos,u,local_grid);

   D        = rho*Lorentz;
   U        = rho*eos.h*Lorentz*Lorentz - p;
   tau      = U - D;
   S_cov[0] = rho*eos.h*Lorentz*Lorentz*v_cov[0];
   S_cov[1] = rho*eos.h*Lorentz*Lorentz*v_cov[1];
   S_cov[2] = rho*eos.h*Lorentz*Lorentz*v_cov[2];

   q[0] = D;
   q[1] = tau;
   q[2] = S_cov[0];
   q[3] = S_cov[1];
   q[4] = S_cov[2];
}
