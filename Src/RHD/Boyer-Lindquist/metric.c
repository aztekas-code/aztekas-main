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

   double Sigma, rho2, Delta, a, M, A;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*cos(theta)*cos(theta);
   Sigma = pow(r*r + a*a,2.0) - a*a*Delta*sin(theta)*sin(theta);

   local_grid->lapse = sqrt(rho2*Delta/Sigma);

   #if POLAR == FALSE

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

   local_grid->dety = sqrt(Sigma*rho2*sin(theta)*sin(theta)/Delta);

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

   double Sigma, rho2, Delta, a, M, A;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*cos(theta)*cos(theta);
   Sigma = pow(r*r + a*a,2.0) - a*a*Delta*sin(theta)*sin(theta);

   #if POLAR == FALSE

   der->dlapse[0] = sqrt(Delta*rho2/Sigma)*(Delta*Sigma*r - Delta*rho2*(a*a*(M - r)*sin(theta) + 2.0*r*(a*a + r*r)) - Sigma*rho2*(M - r))/(Delta*Sigma*rho2);
   der->dlapse[1] = (1.0/2.0)*a*a*sqrt(Delta*rho2/Sigma)*(Delta*rho2 - 2.0*Sigma*sin(theta))*cos(theta)/(Sigma*rho2);
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 2.0*M*a*(a*a*a*a*sin(theta) - a*a*a*a - a*a*r*r*sin(theta) + 2.0*a*a*r*r + 3.0*r*r*r*r)/(Sigma*Sigma);
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = -2.0*Delta*M*a*a*a*r*cos(theta)/(Sigma*Sigma);
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 2.0*(Delta*r + rho2*(M - r))/(Delta*Delta);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*(-Sigma*r + rho2*(a*a*(M - r)*sin(theta) + 2.0*r*(a*a + r*r)))*pow(sin(theta), 2)/(rho2*rho2);

   der->dgam[1][0][0] = -a*a*sin(2.0*theta)/Delta;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = -a*a*sin(2.0*theta);
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = (2.0*Sigma*a*a*pow(sin(theta), 2.0) - rho2*(3.0*Delta*a*a*sin(theta) - 2.0*pow(a*a + r*r, 2.0)))*sin(theta)*cos(theta)/(rho2*rho2);

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

   der->dlapse[0] = sqrt(Delta*rho2/Sigma)*(Delta*Sigma*r - Delta*rho2*(a*a*(M - r)*sin(theta) + 2.0*r*(a*a + r*r)) - Sigma*rho2*(M - r))/(Delta*Sigma*rho2);
   der->dlapse[1] = 0.0;
   der->dlapse[2] = (1.0/2.0)*a*a*sqrt(Delta*rho2/Sigma)*(Delta*rho2 - 2.0*Sigma*sin(theta))*cos(theta)/(Sigma*rho2);

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 2.0*M*a*(a*a*a*a*sin(theta) - a*a*a*a - a*a*r*r*sin(theta) + 2.0*a*a*r*r + 3.0*r*r*r*r)/(Sigma*Sigma);
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = -2.0*Delta*M*a*a*a*r*cos(theta)/(Sigma*Sigma);
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = 2.0*(Delta*r + rho2*(M - r))/(Delta*Delta);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*(-Sigma*r + rho2*(a*a*(M - r)*sin(theta) + 2.0*r*(a*a + r*r)))*pow(sin(theta), 2)/(rho2*rho2);
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*r;

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 0.0;

   der->dgam[2][0][0] = -a*a*sin(2.0*theta)/Delta;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = (2.0*Sigma*a*a*pow(sin(theta), 2.0) - rho2*(3.0*Delta*a*a*sin(theta) - 2.0*pow(a*a + r*r, 2.0)))*sin(theta)*cos(theta)/(rho2*rho2);
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = -a*a*sin(2.0*theta);

   #endif
#endif
}
