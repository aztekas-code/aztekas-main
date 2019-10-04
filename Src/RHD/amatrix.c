#include"main.h"
    
void Matrix_A(double *a, double *u, gauge_ local_grid)
{
   double rho, p, v_cov[3], v_con[3];
   double VV, Lorentz, h;
   double dLorentzdv1, dLorentzdv2, dLorentzdv3;
   double vdLorentz, denominator;
   eos_ eos;

   // Density and Pressure
   rho = u[0];
   p   = u[1];

   // Covariant components of the 3-velocity
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
   Lorentz = 1.0/sqrt(1.0 - VV);

   // Equation of State
   EoS(&eos,u,local_grid);
   h  = eos.h;

   dLorentzdv1 = pow(Lorentz,3.0)*(v_con[0]*v_con[0]);
   dLorentzdv2 = pow(Lorentz,3.0)*(v_con[1]*v_con[1]);
   dLorentzdv3 = pow(Lorentz,3.0)*(v_con[2]*v_con[2]);

   vdLorentz = v_cov[0]*dLorentzdv1 + v_cov[1]*dLorentzdv2 + v_cov[2]*dLorentzdv3;

   denominator = rho*eos.dhdrho*vdLorentz \
               + rho*pow(Lorentz,3.0)*h*eos.dhdp \
               - h*vdLorentz \
               - h*Lorentz;

   a(0,0) = (rho*eos.dhdp*vdLorentz \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - 2.0*h*vdLorentz \
          - h*Lorentz) / (denominator*Lorentz);
   a(0,1) = (rho*eos.dhdp*vdLorentz) / (denominator);
   a(0,2) = (rho*pow(Lorentz,2.0)*eos.dhdp - 1.0)*dLorentzdv1 / (denominator*pow(Lorentz,2.0));

   a(1,0) = (h*vdLorentz \
          - pow(Lorentz*h,2.0) \
          + h*Lorentz \
          - rho*eos.dhdrho*vdLorentz \
          - rho*h*pow(Lorentz,2.0)*eos.dhdrho) / (denominator);
   a(1,1) = (h*vdLorentz \
          + h*Lorentz \
          - rho*eos.dhdrho*vdLorentz) / (denominator);
   a(1,2) = ((rho*eos.dhdrho - h)*dLorentzdv1) / (denominator);

   a(2,0) = v_cov[0]*(rho*eos.dhdrho + h - rho*Lorentz*eos.dhdp) / (rho*denominator);
   a(2,1) = v_cov[0]*Lorentz*eos.dhdp / (denominator);
   a(2,2) = (rho*eos.dhdrho*(vdLorentz - v_cov[0]*dLorentzdv1) \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - h*(vdLorentz - v_cov[0]*dLorentzdv1) \
          - h*Lorentz) / (rho*denominator*pow(Lorentz,2.0)*h);
/*
   a(0,0) = (rho*eos.dhdp*vdLorentz \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - 2.0*h*vdLorentz \
          - h*Lorentz) / (denominator*Lorentz);
   a(0,1) = (rho*eos.dhdp*vdLorentz) / (denominator);
   a(0,2) = (rho*pow(Lorentz,2.0)*eos.dhdp - 1.0)*dLorentzdv1 / (denominator*pow(Lorentz,2.0));
   a(0,3) = (rho*pow(Lorentz,2.0)*eos.dhdp - 1.0)*dLorentzdv2 / (denominator*pow(Lorentz,2.0));
   a(0,4) = (rho*pow(Lorentz,2.0)*eos.dhdp - 1.0)*dLorentzdv3 / (denominator*pow(Lorentz,2.0));

   a(1,0) = (h*vdLorentz \
          - pow(Lorentz*h,2.0) \
          + h*Lorentz \
          - rho*eos.dhdrho*vdLorentz \
          - rho*h*pow(Lorentz,2.0)*eos.dhdrho) / (denominator);
   a(1,1) = (h*vdLorentz \
          + h*Lorentz \
          - rho*eos.dhdrho*vdLorentz) / (denominator);
   a(1,2) = ((rho*eos.dhdrho - h)*dLorentzdv1) / (denominator);
   a(1,3) = ((rho*eos.dhdrho - h)*dLorentzdv1) / (denominator);
   a(1,4) = ((rho*eos.dhdrho - h)*dLorentzdv1) / (denominator);

   a(2,0) = v_cov[0]*(rho*eos.dhdrho + h - rho*Lorentz*eos.dhdp) / (rho*denominator);
   a(2,1) = v_cov[0]*Lorentz*eos.dhdp / (denominator);
   a(2,2) = (rho*eos.dhdrho*(vdLorentz - v_cov[0]*dLorentzdv1) \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - h*(vdLorentz - v_cov[0]*dLorentzdv1) \
          - h*Lorentz) / (rho*denominator*pow(Lorentz,2.0)*h);
   a(2,3) = v_cov[0]*(rho*eos.dhdrho - h)*dLorentzdv2 / (rho*denominator*pow(Lorentz,2.0)*h);
   a(2,4) = v_cov[0]*(rho*eos.dhdrho - h)*dLorentzdv3 / (rho*denominator*pow(Lorentz,2.0)*h);

   a(3,0) = v_cov[1]*(rho*eos.dhdrho + h - rho*Lorentz*eos.dhdp) / (rho*denominator);
   a(3,1) = v_cov[1]*Lorentz*eos.dhdp / (denominator);
   a(3,2) = v_cov[1]*(rho*eos.dhdrho - h)*dLorentzdv1 / (rho*denominator*pow(Lorentz,2.0)*h);
   a(3,3) = (rho*eos.dhdrho*(vdLorentz - v_cov[1]*dLorentzdv2) \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - h*(vdLorentz - v_cov[1]*dLorentzdv2) \
          - h*Lorentz) / (rho*denominator*pow(Lorentz,2.0)*h);
   a(3,4) = v_cov[1]*(rho*eos.dhdrho - h)*dLorentzdv3 / (rho*denominator*pow(Lorentz,2.0)*h);

   a(4,0) = v_cov[2]*(rho*eos.dhdrho + h - rho*Lorentz*eos.dhdp) / (rho*denominator);
   a(4,1) = v_cov[2]*Lorentz*eos.dhdp / (denominator);
   a(4,2) = v_cov[2]*(rho*eos.dhdrho - h)*dLorentzdv1 / (rho*denominator*pow(Lorentz,2.0)*h);
   a(4,3) = v_cov[2]*(rho*eos.dhdrho - h)*dLorentzdv3 / (rho*denominator*pow(Lorentz,2.0)*h);
   a(4,4) = (rho*eos.dhdrho*(vdLorentz - v_cov[2]*dLorentzdv3) \
          + rho*pow(Lorentz,3.0)*h*eos.dhdp \
          - h*(vdLorentz - v_cov[2]*dLorentzdv3) \
          - h*Lorentz) / (rho*denominator*pow(Lorentz,2.0)*h);
*/
}
