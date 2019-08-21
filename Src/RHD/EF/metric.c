#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   local_grid->lapse = sqrt(1.0 - 2.0*MM/r);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0 - 2.0*MM*x*x/(r*r*r);
   local_grid->gamma_con[0][1] = -2.0*MM*x*y/(r*r*r);
   local_grid->gamma_con[0][2] = -2.0*MM*x*z/(r*r*r);
   local_grid->gamma_con[1][0] = -2.0*MM*x*y/(r*r*r);
   local_grid->gamma_con[1][1] = 1.0 - 2.0*MM*y*y/(r*r*r);
   local_grid->gamma_con[1][2] = -2.0*MM*y*z/(r*r*r);
   local_grid->gamma_con[2][0] = -2.0*MM*x*z/(r*r*r);
   local_grid->gamma_con[2][1] = -2.0*MM*y*z/(r*r*r);
   local_grid->gamma_con[2][2] = 1.0 - 2.0*MM*z*z/(r*r*r);

   local_grid->dety = sqrt(r/(r - 2.0*MM));

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];
   double z = local_grid->x[2];
   double r = sqrt(R*R + z*z);

   local_grid->lapse = sqrt(1.0 - 2.0*MM/r);

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1 - 2.0*MM*R*R/(r*r*r);
   local_grid->gamma_con[0][1] = -2.0*MM*R*z/(r*r*r);
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = -2.0*MM*R*z/(r*r*r);
   local_grid->gamma_con[1][1] = 1.0 - 2.0*MM*z*z/(r*r*r);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(R*R);

   local_grid->dety = sqrt(R*R*r/(r - 2.0*MM));

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 

   local_grid->lapse = sqrt(r/(2*MM + r));

   local_grid->beta_con[0] = 2.0*MM/(2.0*MM + r);
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = r/(2.0*MM + r);
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/(r*r);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(r*r*sin(theta)*sin(theta));

   local_grid->dety = pow(r,3.0/2.0)*fabs(sin(theta))*sqrt(2*MM + r);

#endif
}

void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   
   der->dlapse[0] = MM*x/(pow(r,5.0/2.0)*sqrt(r - 2.0*MM));
   der->dlapse[1] = MM*y/(pow(r,5.0/2.0)*sqrt(r - 2.0*MM));
   der->dlapse[2] = MM*z/(pow(r,5.0/2.0)*sqrt(r - 2.0*MM));

   der->dbeta[0][0] = 0.0;
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = -2.0*MM*x*(4.0*MM*r*r - 4.0*MM*x*x - 2.0*r*r*r + 3.0*r*x*x)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][1] = -2.0*MM*y*(r*r*(2.0*MM - r) + r*x*x - 2.0*x*x*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][0][2] = -2.0*MM*z*(r*r*(2.0*MM - r) + r*x*x - 2.0*x*x*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][0] = -2.0*MM*y*(r*r*(2.0*MM - r) + r*x*x - 2.0*x*x*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][1] = 2.0*MM*x*(4.0*MM*y*y - 3.0*r*y*y)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][1][2] = 2.0*MM*x*y*z*(4.0*MM - 3.0*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][0] = -2.0*MM*z*(r*r*(2.0*MM - r) + r*x*x - 2.0*x*x*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][1] = 2.0*MM*x*y*z*(4.0*MM - 3.0*r)/(pow(r,4.0)*pow(2.0*M - r,2.0));
   der->dgam[0][2][2] = 2.0*MM*x*(4.0*MM*z*z - 3.0*r*z*z)/(pow(r,4.0)*pow(2.0*M - r,2.0));

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

   der->dlapse[0] = MM*R/(pow(r,5.0/2.0)*sqrt(r - 2.0*MM));
   der->dlapse[1] = MM*z/(pow(r,5.0/2.0)*sqrt(r - 2.0*MM));
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

   der->dgam[0][0][0] = 2.0*MM*R*(3.0*r*z*z - 4.0*MM*z*z - r*r*r)/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[0][0][1] = -2.0*MM*z*(R*R*r - 2.0*R*R*(2.0*MM - r) + r*r*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[0][0][2] = 0.0;
   der->dgam[0][1][0] = -2.0*MM*z*(R*R*r - 2.0*R*R*(2.0*MM - r) + r*r*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[0][1][1] = -2.0*MM*R(4.0*MM*R*R - 4.0*M*r*r - 3.0*R*R*r + 3.0*r*r*r)/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 0.0;
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*R;

   der->dgam[1][0][0] = 2.0*MM*z*(4.0*MM*r*r - 4.0*MM*z*z - 3.0*r*r*r + 3.0*r*z*z)/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[1][0][1] = -2.0*MM*R*(r*r*(2.0*MM - r) + r*z*z - 2.0*z*z*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[1][0][2] = 0.0;
   der->dgam[1][1][0] = -2.0*MM*R*(r*r*(2.0*MM - r) + r*z*z - 2.0*z*z*(2.0*MM - r))/(pow(r,4.0)*pow(2.0*MM - r,2.0));
   der->dgam[1][1][1] = 2.0*MM*z*(3.0*r*R*R - 4.0*MM*R*R - r*r*r)/(pow(r,4.0)*pow(2.0*MM - r,2.0));
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

   der->dlapse[0] = MM/(sqrt(r*pow(2.0*MM + r,3.0)));
   der->dlapse[1] = 0.0;
   der->dlapse[2] = 0.0;

   der->dbeta[0][0] = -2.0*MM/pow(2.0*MM + r,2.0);
   der->dbeta[0][1] = 0.0;
   der->dbeta[0][2] = 0.0;
   der->dbeta[1][0] = 0.0;
   der->dbeta[1][1] = 0.0;
   der->dbeta[1][2] = 0.0;
   der->dbeta[2][0] = 0.0;
   der->dbeta[2][1] = 0.0;
   der->dbeta[2][2] = 0.0;

   der->dgam[0][0][0] = -2.0*MM/(r*r);
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
