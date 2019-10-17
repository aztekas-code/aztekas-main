#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   double M = Black_Hole_Mass;
   
   local_grid->lapse = sqrt(1.0 - 2.0*M/r);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0 - 2.0*M*x*x/(r*r*r);
   local_grid->gamma_con[0][1] = -2.0*M*x*y/(r*r*r);
   local_grid->gamma_con[0][2] = -2.0*M*x*z/(r*r*r);
   local_grid->gamma_con[1][0] = -2.0*M*x*y/(r*r*r);
   local_grid->gamma_con[1][1] = 1.0 - 2.0*M*y*y/(r*r*r);
   local_grid->gamma_con[1][2] = -2.0*M*y*z/(r*r*r);
   local_grid->gamma_con[2][0] = -2.0*M*x*z/(r*r*r);
   local_grid->gamma_con[2][1] = -2.0*M*y*z/(r*r*r);
   local_grid->gamma_con[2][2] = 1.0 - 2.0*M*z*z/(r*r*r);

   local_grid->dety = sqrt(r/(r - 2.0*M));

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);
   double M = Black_Hole_Mass;

   local_grid->lapse = sqrt(1.0 - 2.0*M/r);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1 - 2.0*M*R*R/(r*r*r);
   local_grid->gamma_con[0][1] = -2.0*M*R*z/(r*r*r);
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = -2.0*M*R*z/(r*r*r);
   local_grid->gamma_con[1][1] = 1.0 - 2.0*M*z*z/(r*r*r);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(R*R);

   local_grid->dety = sqrt(R*R*r/(r - 2.0*M));

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 
   double M = Black_Hole_Mass;

   local_grid->lapse = sqrt(1.0 - 2.0*M/r);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0 - 2.0*M/r;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/(r*r);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(r*r*sin(theta)*sin(theta));

   local_grid->dety = pow(r,5.0/2.0)*fabs(sin(theta))*sqrt(1.0/(r - 2.0*M));

#endif
}

void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   double M = Black_Hole_Mass;
   
   der->dlapse[0] = M*x/(pow(r,5.0/2.0)*sqrt(r - 2.0*M));
   der->dlapse[1] = M*y/(pow(r,5.0/2.0)*sqrt(r - 2.0*M));
   der->dlapse[2] = M*z/(pow(r,5.0/2.0)*sqrt(r - 2.0*M));

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = -2.0*M*x*(4.0*M*r*r - 4.0*M*x*x - 2.0*r*r*r + 3.0*r*x*x)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][1] = -2.0*M*y*(r*r*(2.0*M - r) + r*x*x - 2.0*x*x*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][2] = -2.0*M*z*(r*r*(2.0*M - r) + r*x*x - 2.0*x*x*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][0] = -2.0*M*y*(r*r*(2.0*M - r) + r*x*x - 2.0*x*x*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][1] = 2.0*M*x*(4.0*M*y*y - 3.0*r*y*y)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][2] = 2.0*M*x*y*z*(4.0*M - 3.0*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][0] = -2.0*M*z*(r*r*(2.0*M - r) + r*x*x - 2.0*x*x*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][1] = 2.0*M*x*y*z*(4.0*M - 3.0*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][2] = 2.0*M*x*(4.0*M*z*z - 3.0*r*z*z)/(pow(r,4.0)*pow(2.0*M - r,2.0));

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
   double M = Black_Hole_Mass;

   der->dlapse[0] = M*R/(pow(r,5.0/2.0)*sqrt(r - 2.0*M));
   der->dlapse[1] = M*z/(pow(r,5.0/2.0)*sqrt(r - 2.0*M));
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

   der->dgam[0][0][0] = 2.0*M*R*(3.0*r*z*z - 4.0*M*z*z - r*r*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][1] = -2.0*M*z*(R*R*r - 2.0*R*R*(2.0*M - r) + r*r*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = -2.0*M*z*(R*R*r - 2.0*R*R*(2.0*M - r) + r*r*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][1] = -2.0*M*R(4.0*M*R*R - 4.0*M*r*r - 3.0*R*R*r + 3.0*r*r*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*R;

   der->dgam[1][0][0] = 2.0*M*z*(4.0*M*r*r - 4.0*M*z*z - 3.0*r*r*r + 3.0*r*z*z)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[1][0][1] = -2.0*M*R*(r*r*(2.0*M - r) + r*z*z - 2.0*z*z*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = -2.0*M*R*(r*r*(2.0*M - r) + r*z*z - 2.0*z*z*(2.0*M - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[1][1][1] = 2.0*M*z*(3.0*r*R*R - 4.0*M*R*R - r*r*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
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
   double M = Black_Hole_Mass;

   der->dlapse[0] = M/(pow(r,3.0/2.0)*sqrt(r - 2.0*M));
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

   der->dgam[0][0][0] = -2.0*M/pow(2.0*M - r,2.0);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*r*sin(theta)*sin(theta);

   der->dgam[1][0][0] = 0.0;
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = 0.0;
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = 0.0;
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = 2.0*r*r*cos(theta)*sin(theta);

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
