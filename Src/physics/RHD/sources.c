/**
 * @file /RHD/sources.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 */

#include"main.h"
    
void Source_Terms(double *s, double *u, gauge_ *local_grid)
{
   int i;
   double rho, p, v_cov[3], v_con[3];
   double D, U, tau, Lorentz, h, cs, VV;
   double S_cov[3], S_con[3], V[3];
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
   v_con[0] = local_grid->gamma_con[0][0]*v_cov[0] + \
              local_grid->gamma_con[0][1]*v_cov[1] + \
              local_grid->gamma_con[0][2]*v_cov[2];
   v_con[1] = local_grid->gamma_con[1][0]*v_cov[0] + \
              local_grid->gamma_con[1][1]*v_cov[1] + \
              local_grid->gamma_con[1][2]*v_cov[2];
   v_con[2] = local_grid->gamma_con[2][0]*v_cov[0] + \
              local_grid->gamma_con[2][1]*v_cov[1] + \
              local_grid->gamma_con[2][2]*v_cov[2];

   // Contraction v_i v^i
   VV = v_con[0]*v_cov[0] + v_con[1]*v_cov[1] + v_con[2]*v_cov[2];

   if(VV > 1.0)
   {
      VV = 0.99;
   }

   // Lorentz Factor
   Lorentz = 1.0/sqrt(1.0 - VV);

   // Equation of State
   EoS(&eos,u,local_grid);
   h  = eos.h;
   cs = eos.cs;

   // Define conservative density D, variable U and conservative energy \tau
   D   = rho*Lorentz;
   U   = rho*h*Lorentz*Lorentz - p;
   tau = U - D;

   // Compute the covariant and contravariant components of the 3-momentum
   for(i = 0; i < 3; i++)
   {
      S_cov[i] = rho*eos.h*Lorentz*Lorentz*v_cov[i];
      S_con[i] = rho*eos.h*Lorentz*Lorentz*v_con[i];
   }

   double W[3][3];
   double Wd[3][3];
   double gam[3][3];
   double beta[3];
   double lapse;
   double Sd[3];
   double Su[3];
   double Wik_dgamik[3];
   double Wik_betaj_dgamjik;
   double lapse_Wik_dgam1ik;
   double lapse_Wik_dgam2ik;
   double lapse_Wik_dgam3ik;
   double Wdij_dbetaji;
   double Sdi_dbeta1i;
   double Sdi_dbeta2i;
   double Sdi_dbeta3i;
   double Suj_dlapsej;
   double U_dlapse1;
   double U_dlapse2;
   double U_dlapse3;
   der_gauge_ der;

   Gauge_Derivatives(&der,local_grid);

   lapse     = local_grid->lapse;
   beta[0]   = local_grid->beta_con[0];
   beta[1]   = local_grid->beta_con[1];
   beta[2]   = local_grid->beta_con[2];
   gam[0][0] = local_grid->gamma_con[0][0];
   gam[0][1] = local_grid->gamma_con[0][1];
   gam[0][2] = local_grid->gamma_con[0][2];
   gam[1][0] = local_grid->gamma_con[1][0];
   gam[1][1] = local_grid->gamma_con[1][1];
   gam[1][2] = local_grid->gamma_con[1][2];
   gam[2][0] = local_grid->gamma_con[2][0];
   gam[2][1] = local_grid->gamma_con[2][1];
   gam[2][2] = local_grid->gamma_con[2][2];

   Su[0] = S_con[0];
   Su[1] = S_con[1];
   Su[2] = S_con[2];

   Sd[0] = S_cov[0];
   Sd[1] = S_cov[1];
   Sd[2] = S_cov[2];

   W[0][0]  = Su[0]*v_con[0] + p*gam[0][0]; 
   W[0][1]  = Su[0]*v_con[1] + p*gam[0][1]; 
   W[0][2]  = Su[0]*v_con[2] + p*gam[0][2]; 
   W[1][0]  = Su[1]*v_con[0] + p*gam[1][0]; 
   W[1][1]  = Su[1]*v_con[1] + p*gam[1][1]; 
   W[1][2]  = Su[1]*v_con[2] + p*gam[1][2]; 
   W[2][0]  = Su[2]*v_con[0] + p*gam[2][0]; 
   W[2][1]  = Su[2]*v_con[1] + p*gam[2][1]; 
   W[2][2]  = Su[2]*v_con[2] + p*gam[2][2]; 

   Wd[0][0] = Sd[0]*v_con[0] + p;
   Wd[0][1] = Sd[0]*v_con[1];
   Wd[0][2] = Sd[0]*v_con[2];
   Wd[1][0] = Sd[1]*v_con[0];
   Wd[1][1] = Sd[1]*v_con[1] + p;
   Wd[1][2] = Sd[1]*v_con[2];
   Wd[2][0] = Sd[2]*v_con[0];
   Wd[2][1] = Sd[2]*v_con[1];
   Wd[2][2] = Sd[2]*v_con[2] + p;

   Wik_dgamik[0] = W[0][0]*der.dgam[0][0][0] + W[0][1]*der.dgam[0][0][1] + W[0][2]*der.dgam[0][0][2] + \
                   W[1][0]*der.dgam[0][1][0] + W[1][1]*der.dgam[0][1][1] + W[1][2]*der.dgam[0][1][2] + \
                   W[2][0]*der.dgam[0][2][0] + W[2][1]*der.dgam[0][2][1] + W[2][2]*der.dgam[0][2][2];

   Wik_dgamik[1] = W[0][0]*der.dgam[1][0][0] + W[0][1]*der.dgam[1][0][1] + W[0][2]*der.dgam[1][0][2] + \
                   W[1][0]*der.dgam[1][1][0] + W[1][1]*der.dgam[1][1][1] + W[1][2]*der.dgam[1][1][2] + \
                   W[2][0]*der.dgam[1][2][0] + W[2][1]*der.dgam[1][2][1] + W[2][2]*der.dgam[1][2][2];

   Wik_dgamik[2] = W[0][0]*der.dgam[2][0][0] + W[0][1]*der.dgam[2][0][1] + W[0][2]*der.dgam[2][0][2] + \
                   W[1][0]*der.dgam[2][1][0] + W[1][1]*der.dgam[2][1][1] + W[1][2]*der.dgam[2][1][2] + \
                   W[2][0]*der.dgam[2][2][0] + W[2][1]*der.dgam[2][2][1] + W[2][2]*der.dgam[2][2][2];

   Wik_betaj_dgamjik = Wik_dgamik[0]*beta[0] + Wik_dgamik[1]*beta[1] + Wik_dgamik[2]*beta[2];

   lapse_Wik_dgam1ik = lapse*Wik_dgamik[0];
   lapse_Wik_dgam2ik = lapse*Wik_dgamik[1];
   lapse_Wik_dgam3ik = lapse*Wik_dgamik[2];

   Wdij_dbetaji = Wd[0][0]*der.dbeta[0][0] + Wd[0][1]*der.dbeta[0][1] + Wd[0][2]*der.dbeta[0][2] + \
                  Wd[1][0]*der.dbeta[1][0] + Wd[1][1]*der.dbeta[1][1] + Wd[1][2]*der.dbeta[1][2] + \
                  Wd[2][0]*der.dbeta[2][0] + Wd[2][1]*der.dbeta[2][1] + Wd[2][2]*der.dbeta[2][2];

   Sdi_dbeta1i = Sd[0]*der.dbeta[0][0] + Sd[1]*der.dbeta[0][1] + Sd[2]*der.dbeta[0][2];
   Sdi_dbeta2i = Sd[0]*der.dbeta[1][0] + Sd[1]*der.dbeta[1][1] + Sd[2]*der.dbeta[1][2];
   Sdi_dbeta3i = Sd[0]*der.dbeta[2][0] + Sd[1]*der.dbeta[2][1] + Sd[2]*der.dbeta[2][2];

   Suj_dlapsej = Su[0]*der.dlapse[0] + Su[1]*der.dlapse[1] + Su[2]*der.dlapse[2];

   U_dlapse1 = U*der.dlapse[0];
   U_dlapse2 = U*der.dlapse[1];
   U_dlapse3 = U*der.dlapse[2];

   s[0] = 0.0;
   s[1] = Wik_betaj_dgamjik/2.0 + Wdij_dbetaji - Suj_dlapsej;
   s[2] = lapse_Wik_dgam1ik/2.0 + Sdi_dbeta1i - U_dlapse1;
   s[3] = lapse_Wik_dgam2ik/2.0 + Sdi_dbeta2i - U_dlapse2;
   s[4] = lapse_Wik_dgam3ik/2.0 + Sdi_dbeta3i - U_dlapse3;
}
