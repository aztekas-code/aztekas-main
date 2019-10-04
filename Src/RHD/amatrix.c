#include"main.h"
    
void Matrix_A(double *a, double *u, gauge_ local_grid)
{
   double rho, p, v_cov[3], v_con[3];
   double VV, L, h;
   double v0, v1, v2;
   double dLdv0, dLdv1, dLdv2;
   double dhdrho, dhdp;
   double vdL, denominator;
   eos_ eos;

   // Density and Pressure
   rho = u[RHO];
   p   = u[PRE];

   // Covariant components of the 3-velocity
#if DIM == 1
   v_cov[0] = u[VX1];
   v_cov[1] = 0.0;
   v_cov[2] = 0.0;
#elif DIM == 2
   v_cov[0] = u[VX1];
   v_cov[1] = u[VX2];
   v_cov[2] = 0.0;
#elif DIM == 3 || DIM == 4
   v_cov[0] = u[VX1];
   v_cov[1] = u[VX2];
   v_cov[2] = u[VX3];
#endif

   // Contravariant components of the 3-velocity
   v_con[0] = local_grid.gamma_con[0][0]*v_cov[0] + \
              local_grid.gamma_con[0][1]*v_cov[1] + \
              local_grid.gamma_con[0][2]*v_cov[2];
   v_con[1] = local_grid.gamma_con[1][0]*v_cov[0] + \
              local_grid.gamma_con[1][1]*v_cov[1] + \
              local_grid.gamma_con[1][2]*v_cov[2];
   v_con[2] = local_grid.gamma_con[2][0]*v_cov[0] + \
              local_grid.gamma_con[2][1]*v_cov[1] + \
              local_grid.gamma_con[2][2]*v_cov[2];

   // Contraction v_i v^i
   VV = v_con[0]*v_cov[0] + v_con[1]*v_cov[1] + v_con[2]*v_cov[2];

   // Lorentz Factor
   L = 1.0/sqrt(1.0 - VV);

   // Equation of State
   EoS(&eos,u,local_grid);
   h      = eos.h;
   dhdrho = eos.dhdrho;
   dhdp   = eos.dhdp;

   v0 = v_cov[0];
   v1 = v_cov[1];
   v2 = v_cov[2];

   dLdv0 = pow(L,3.0)*(v_con[0]);
   dLdv1 = pow(L,3.0)*(v_con[1]);
   dLdv2 = pow(L,3.0)*(v_con[2]);

   vdL = v_cov[0]*dLdv0 + v_cov[1]*dLdv1 + v_cov[2]*dLdv2;

   denominator = rho*dhdrho*vdL \
               + rho*pow(L,3.0)*h*dhdp \
               - h*vdL \
               - h*L;

   a(0,0) = (rho*dhdp*L*vdL \
          + rho*pow(L,3.0)*h*dhdp \
          - 2.0*h*vdL \
          - h*L) / (denominator*L);
   a(0,1) = (rho*dhdp*vdL) / (denominator);
   a(0,2) = -(rho*pow(L,2.0)*dhdp - 1.0)*dLdv0 / (denominator*pow(L,2.0));
   a(0,3) = -(rho*pow(L,2.0)*dhdp - 1.0)*dLdv1 / (denominator*pow(L,2.0));
   a(0,4) = -(rho*pow(L,2.0)*dhdp - 1.0)*dLdv2 / (denominator*pow(L,2.0));

   a(1,0) = (-rho*dhdrho*vdL \
          - rho*pow(L,2.0)*h*dhdrho \
          + h*vdL \
          - pow(L*h,2.0) \
          + L*h) / (denominator);
   a(1,1) = (-rho*dhdrho*vdL \
          + h*vdL \
          + L*h) / (denominator);
   a(1,2) = ((rho*dhdrho - h)*dLdv0) / (denominator);
   a(1,3) = ((rho*dhdrho - h)*dLdv0) / (denominator);
   a(1,4) = ((rho*dhdrho - h)*dLdv0) / (denominator);

   a(2,0) = v_cov[0]*(rho*dhdrho + h - rho*L*dhdp) / (rho*denominator);
   a(2,1) = -v_cov[0]*L*dhdp / (denominator);
   a(2,2) = (rho*dhdrho*(vdL - v_cov[0]*dLdv0) \
          + rho*pow(L,3.0)*h*dhdp \
          - h*(vdL - v_cov[0]*dLdv0) \
          - h*L) / (rho*denominator*pow(L,2.0)*h);
   a(2,3) = -v_cov[0]*(rho*dhdrho - h)*dLdv1 / (rho*denominator*pow(L,2.0)*h);
   a(2,4) = -v_cov[0]*(rho*dhdrho - h)*dLdv2 / (rho*denominator*pow(L,2.0)*h);

   a(3,0) = v_cov[1]*(rho*dhdrho + h - rho*L*dhdp) / (rho*denominator);
   a(3,1) = -v_cov[1]*L*dhdp / (denominator);
   a(3,2) = -v_cov[1]*(rho*dhdrho - h)*dLdv0 / (rho*denominator*pow(L,2.0)*h);
   a(3,3) = (rho*dhdrho*(vdL - v_cov[1]*dLdv1) \
          + rho*pow(L,3.0)*h*dhdp \
          - h*(vdL - v_cov[1]*dLdv1) \
          - h*L) / (rho*denominator*pow(L,2.0)*h);
   a(3,4) = -v_cov[1]*(rho*dhdrho - h)*dLdv2 / (rho*denominator*pow(L,2.0)*h);

   a(4,0) = v_cov[2]*(rho*dhdrho + h - rho*L*dhdp) / (rho*denominator);
   a(4,1) = -v_cov[2]*L*dhdp / (denominator);
   a(4,2) = -v_cov[2]*(rho*dhdrho - h)*dLdv0 / (rho*denominator*pow(L,2.0)*h);
   a(4,3) = -v_cov[2]*(rho*dhdrho - h)*dLdv2 / (rho*denominator*pow(L,2.0)*h);
   a(4,4) = (rho*dhdrho*(vdL - v_cov[2]*dLdv2) \
          + rho*pow(L,3.0)*h*dhdp \
          - h*(vdL - v_cov[2]*dLdv2) \
          - h*L) / (rho*denominator*pow(L,2.0)*h);

}
