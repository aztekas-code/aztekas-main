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

   double Sigma, Delta, a, M, A;
   double term1, term2;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*pow(cos(theta),2.0);
   Sigma = por(r*r + a*a,2.0) - a*a*Delta*pow(sin(theta),2.0);

   local_grid->lapse = 1.0/sqrt(1.0 + 2.0*M*r/rho2);

   local_grid->beta_con[0] = (2.0*M*r/rho2)/(1.0 + 2.0*M*r/rho2)
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   #if POLAR == FALSE

   local_grid->gamma_con[0][0] = (a*a*(rho2 + 2.0*M*r)*pow(sin(theta),2.0) + rho2*rho2)/(rho2*(rho2 + 2.0*M*r));
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = a/rho2;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/rho2;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = a/rho2;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(rho2*pow(sin(theta),2.0));

   #elif POLAR == TRUE

   local_grid->gamma_con[0][0] = (a*a*(rho2 + 2.0*M*r)*pow(sin(theta),2.0) + rho2*rho2)/(rho2*(rho2 + 2.0*M*r));
   local_grid->gamma_con[0][1] = a/rho2;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = a/rho2;
   local_grid->gamma_con[1][1] = 1.0/(rho2*pow(sin(theta),2.0));
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/rho2;

   #endif

   local_grid->dety = sqrt(local_grid->lapse*rho2*(2.0*M*r + rho2)*pow(sin(theta),2.0));

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

   double Sigma, Delta, a, M, A;
   double term1, term2, term3;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Delta = r*r - 2.0*M*r + a*a;
   rho2  = r*r + a*a*pow(cos(theta),2.0);
   Sigma = por(r*r + a*a,2.0) - a*a*Delta*pow(sin(theta),2.0);

   der->dlapse[0] = (M*(r*r - a*a*pow(cos(theta),2.0)))/(rho2*rho2*pow(1.0 + 2.0*M*r/rho2,1.5));
   der->dlapse[1] = -(2.0*M*a*a*r*sin(theta)*cos(theta))/(rho2*rho2*pow(1.0 + 2.0*M*r/rho2,1.5)); 
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = (2.0*(a*a*pow(cos(theta),2.0) - r*r)/pow(rho2 + 2.0*M*r,2.0);
   der->dbeta[0][1] = 8*M*a*a*r*sin(2.0*theta)/pow(2.0*r*r + a*a + a*a*cos(2.0*theta) + 4.0*M*r);
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   #if POLAR == FALSE

   der->dgam[0][0][0] = 2.0*M*(a*a*pow(cos(theta),2.0 - r*r)/pow(rho2,2.0);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 2.0*M*a*(r*r - a*a*pow(cos(theta),2.0))*pow(sin(theta),2.0)/pow(rho2,2.0);
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 2.0*M*a*(r*r - a*a*pow(cos(theta),2.0))*pow(sin(theta),2.0)/pow(rho2,2.0);
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*(r*rho2*rho2 - M*a*a*(r*r - a*a*pow(cos(theta),2.0))*pow(sin(theta),2.0))*pow(sin(theta),2.0)/pow(rho2,2.0);

   der->dgam[1][0][0] = 4.0*M*a*a*r*sin(theta)*cos(theta)/pow(rho2,2.0);
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = -2.0*a*(2.0*M*a*a*r*pow(sin(theta),2.0) + rho2*(2.0*M*r + rho2))*sin(theta)*cos(theta)/pow(rho2,2.0);
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = -2.0*a*a*sin(theta)*cos(theta);
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = -2.0*a*(2.0*M*a*a*r*pow(sin(theta),2.0) + rho2*(2.0*M*r + rho2))*sin(theta)*cos(theta)/pow(rho2,2.0);
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 2.0*(2.0*M*a*a*r*(a*a + r*r)*pow(sin(theta),2.0) + rho2*(a*a*(rho2 + 2.0*M*r)*pow(sin(theta),2.0) + rho2*rho2))*sin(theta)*cos(theta)/pow(rho2,2.0);

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

   der->dgam[0][0][0] = 2.0*M*(Sigma - 2.0*r*r)/(Sigma*Sigma);
   der->dgam[0][0][1] = 2.0*a*M*(2.0*r*r - Sigma)*sin(theta)*sin(theta)/(Sigma*Sigma);
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 2.0*a*M*(2.0*r*r - Sigma)*sin(theta)*sin(theta)/(Sigma*Sigma);
   der->dgam[0][1][1] = 2.0*M*a*a*pow(sin(theta),4.0)/Sigma - 4.0*M*a*a*r*r*pow(sin(theta),4.0)/(Sigma*Sigma) + 2.0*r*sin(theta)*sin(theta);
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*r;

   #endif

#endif
}
