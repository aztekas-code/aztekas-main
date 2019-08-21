#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN
   
   local_grid->lapse = 1.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0;

   local_grid->dety = 1.0;

#elif COORDINATES == CYLINDRICAL

   double R = local_grid->x[1];

   local_grid->lapse = 1.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(R*R);

   local_grid->dety = R;

#elif COORDINATES == SPHERICAL

   double r     = local_grid->x[1];
   double theta = local_grid->x[2]; 

   local_grid->lapse = 1.0;

   local_grid->beta_con[0] = 0.0;
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = 1.0;
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = 0.0;
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/(r*r);
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = 0.0;
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = 1.0/(r*r*sin(theta)*sin(theta));

   local_grid->dety = r*r*sin(theta);

#endif
}

void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN
   
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
   der->dgam[0][2][2] = 2.0*R;

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
