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

   double rho2, Delta, a, M, A;
   double term1, term2;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   rho2  = r*r + a*a*cos(theta)*cos(theta);
   Delta = r*r - 2.0*M*r + a*a;
   A     = rho2*rho2 + a*a*sin(theta)*sin(theta)*(rho2 + 2.0*M*r);
   term1 = 2.0*M*r + rho2;
   term2 = A - a*a*term1*sin(theta)*sin(theta);

   local_grid->lapse = 1.0/sqrt(1.0 + 2.0*M*r/rho2);

   local_grid->beta_con[0] = (2.0*M*r/rho2)/(1.0 + 2*M*r/rho2);
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   #if POLAR == FALSE

   local_grid->gamma_con[0][0] = A*rho2/(term2*term1);
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = rho2*a/(term2);
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/rho2;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = rho2*a/(term2);
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = rho2/(term2*sin(theta)*sin(theta));

   #elif POLAR == TRUE

   local_grid->gamma_con[0][0] = A*rho2/(term2*term1);
   local_grid->gamma_con[0][1] = rho2*a/(term2);
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = rho2*a/(term2);
   local_grid->gamma_con[1][1] = rho2/(term2*sin(theta)*sin(theta));
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/rho2;

   #endif

   //local_grid->dety = sqrt(term2*term1*sin(theta)*sin(theta)/rho2);
   local_grid->dety = rho2*sin(theta);

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

   double rho2, Delta, a, M, A;
   double term1, term2, term3;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   rho2 = r*r + a*a*cos(theta)*cos(theta);
   Delta = r*r - 2.0*M*r + a*a;
   A     = rho2*rho2 + a*a*sin(theta)*sin(theta)*(rho2 + 2.0*M*r);
   term1 = 2.0*M*r + rho2;
   term2 = A - a*a*term1*sin(theta)*sin(theta);
   term3 = 1 + 2.0*M*r/rho2;

   der->dlapse[0] = (r*r - a*a*cos(theta)*cos(theta))*M/(sqrt(rho2)*pow(term1,3.0/2.0));
   der->dlapse[1] = -M*a*a*r*sin(2.0*theta)/(rho2*rho2*pow(term3,3.0/2.0));
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = 2.0*M*(a*a*cos(theta)*cos(theta) - r*r)/(term1*term1);
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 2.0*M*a*a*r*sin(2.0*theta)/(term1*term1);
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   #if POLAR == FALSE

   der->dgam[0][0][0] = 2.0*M*(rho2 - 2.0*r*r)/(rho2*rho2);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 2.0*a*M*(2.0*r*r - rho2)*sin(theta)*sin(theta)/(rho2*rho2);
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 2.0*a*M*(2.0*r*r - rho2)*sin(theta)*sin(theta)/(rho2*rho2);
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*M*a*a*pow(sin(theta),4.0)/rho2 - 4.0*M*a*a*r*r*pow(sin(theta),4.0)/(rho2*rho2) + 2.0*r*sin(theta)*sin(theta);

   #elif POLAR == TRUE

   der->dgam[0][0][0] = 2.0*M*(rho2 - 2.0*r*r)/(rho2*rho2);
   der->dgam[0][0][1] = 2.0*a*M*(2.0*r*r - rho2)*sin(theta)*sin(theta)/(rho2*rho2);
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 2.0*a*M*(2.0*r*r - rho2)*sin(theta)*sin(theta)/(rho2*rho2);
   der->dgam[0][1][1] = 2.0*M*a*a*pow(sin(theta),4.0)/rho2 - 4.0*M*a*a*r*r*pow(sin(theta),4.0)/(rho2*rho2) + 2.0*r*sin(theta)*sin(theta);
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*r;

   #endif

   der->dgam[1][0][0] = 2.0*M*a*a*r*sin(2.0*theta)/(rho2*rho2);
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = -a*(2.0*M*rho2*r + 2.0*M*a*a*r*sin(theta)*sin(theta) + rho2*rho2)*sin(2.0*theta);
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = -a*a*sin(2.0*theta);
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = -a*(2.0*M*rho2*r + 2.0*M*a*a*r*sin(theta)*sin(theta) + rho2*rho2)*sin(2.0*theta);
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = (4.0*M*rho2*a*a*r*sin(theta)*sin(theta) + 2.0*M*a*a*a*a*r*pow(sin(theta),4.0) + pow(rho2,3.0) + rho2*rho2*a*a*sin(theta)*sin(theta))*sin(2.0*theta)/(rho2*rho2);

   der->dgam[2][0][0] = 0.0;
   der->dgam[2][0][1] = 0.0;
   der->dgam[2][0][2] = 0.0;
   der->dgam[2][1][0] = 0.0;
   der->dgam[2][1][1] = 0.0;
   der->dgam[2][1][2] = 0.0;
   der->dgam[2][2][0] = 0.0;
   der->dgam[2][2][1] = 0.0;
   der->dgam[2][2][2] = 0.0;

#endif
}
