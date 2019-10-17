#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   local_grid->lapse = 0.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 0.0; 
   local_grid->gamma_con[0][1] = 0.0; 
   local_grid->gamma_con[0][2] = 0.0; 
   local_grid->gamma_con[1][0] = 0.0; 
   local_grid->gamma_con[1][1] = 0.0; 
   local_grid->gamma_con[1][2] = 0.0; 
   local_grid->gamma_con[2][0] = 0.0; 
   local_grid->gamma_con[2][1] = 0.0; 
   local_grid->gamma_con[2][2] = 0.0; 

   local_grid->dety = 0.0;

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);

   local_grid->lapse = 0.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 0.0;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 0.0;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 0.0;

   local_grid->dety = 0.0;

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 
   double phi   = local_grid->x[3];

   double Sigma, rho2, Delta, a, M;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*cos(theta)*cos(theta);
   Sigma = pow(r*r + a*a,2.0) - a*a*Delta*sin(theta)*sin(theta);

   #if POLAR == FALSE

   local_grid->lapse = sqrt(rho2*Delta/Sigma);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = -2.0*M*a*r/Sigma;

   local_grid->gamma_con[0][0] = Delta/rho2;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/rho2;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = rho2/(Sigma*sin(theta)*sin(theta));

   #elif POLAR == TRUE

   local_grid->lapse = sqrt(rho2*Delta/Sigma);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = -2*M*a*r/Sigma;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = Delta/rho2;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = rho2/(Sigma*sin(theta)*sin(theta));
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/rho2;

   #endif

   local_grid->dety = sqrt(Sigma*rho2*rho2*sin(theta)*sin(theta)/Delta);

#endif
}

void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   der->dlapse[0] = 0.0;
   der->dlapse[1] = 0.0;
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 0.0;
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 0.0;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 0.0;

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);

   der->dlapse[0] = 0.0;
   der->dlapse[1] = 0.0;
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 0.0;
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 0.0;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 0.0;

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 
   double phi   = local_grid->x[3];

   double Sigma, rho2, Delta, a, M;
   double lapse;
   double drho2dr, drho2dt;
   double dSigmadr, dSigmadt;
   double dDeltadr;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*cos(theta)*cos(theta);
   Sigma = pow(r*r + a*a,2.0) - a*a*Delta*sin(theta)*sin(theta);

   lapse    =  sqrt(rho2*Delta/Sigma);
   drho2dr  =  2.0*r;
   drho2dt  = -2.0*a*a*sin(theta)*cos(theta);
   dSigmadr =  4.0*r*(r*r + a*a) - 2.0*a*a*(r-M)*sin(theta)*sin(theta);
   dSigmadt = -2.0*a*a*Delta*sin(theta)*cos(theta);
   dDeltadr =  2.0*(r-M);

   #if POLAR == FALSE

   der->dlapse[0] = (1.0/(2*lapse))*((drho2dr*Delta*Sigma + rho2*dDeltadr*Sigma - rho2*Delta*dSigmadr)/(Sigma*Sigma));
   der->dlapse[1] = (1.0/(2*lapse))*((drho2dt*Delta*Sigma - rho2*Delta*dSigmadt)/(Sigma*Sigma));
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = -2.0*M*a*((Sigma - r*dSigmadr)/(Sigma*Sigma));
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 2.0*M*a*r*dSigmadt/(Sigma*Sigma);
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = (drho2dr*Delta - rho2*dDeltadr)/(Delta*Delta);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = drho2dr;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = sin(theta)*sin(theta)*(dSigmadr*rho2 - Sigma*drho2dr)/(rho2*rho2);

   der->dgam[1][0][0] = drho2dt/Delta;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = drho2dt;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 2.0*sin(theta)*cos(theta)*(Sigma/rho2) + sin(theta)*sin(theta)*(dSigmadt*rho2 - Sigma*drho2dt)/(rho2*rho2);

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

   #elif POLAR == TRUE

   der->dlapse[0] = (1.0/(2*lapse))*((drho2dr*Delta*Sigma + rho2*dDeltadr*Sigma - rho2*Delta*dSigmadr)/(Sigma*Sigma));
   der->dlapse[2] = (1.0/(2*lapse))*((drho2dt*Delta*Sigma - rho2*Delta*dSigmadt)/(Sigma*Sigma));
   der->dlapse[1] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[0][1] = -2.0*M*a*((Sigma - r*dSigmadr)/(Sigma*Sigma));
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][2] = 0.0;
   der->dbeta[2][1] = 2.0*M*a*r*dSigmadt/(Sigma*Sigma);
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[1][1] = 0.0;

   der->dgam[0][0][0] = (drho2dr*Delta - rho2*dDeltadr)/(Delta*Delta);
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][2] = drho2dr;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][1][1] = sin(theta)*sin(theta)*(dSigmadr*rho2 - Sigma*drho2dr)/(rho2*rho2);

   der->dgam[2][0][0] = (drho2dt*Delta)/(Delta*Delta);
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][2] = drho2dt;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][1][1] = 2.0*sin(theta)*cos(theta)*(Sigma/rho2) + sin(theta)*sin(theta)*(dSigmadr*rho2 - Sigma*drho2dr)/(rho2*rho2);

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   #endif 
   
#endif
}
