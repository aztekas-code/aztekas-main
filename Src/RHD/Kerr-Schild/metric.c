#include"main.h"

void Get_Metric_Components(gauge_ *local_grid)
{
#if COORDINATES == CARTESIAN

   double x = local_grid->x[1];
   double y = local_grid->x[2];
   double z = local_grid->x[3];
   double r = sqrt(x*x + y*y + z*z);
   double M = Black_Hole_Mass;
   
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
   double M = Black_Hole_Mass;

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
   double phi   = local_grid->x[3];
   double Sigma, Delta, a, M, A;
   double term1, term2;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Sigma = r*r + a*a*cos(theta)*cos(theta);
   Delta = r*r - 2.0*M*r + a*a;
   A     = Sigma*Sigma + a*a*sin(theta)*sin(theta)*(Sigma + 2.0*M*r);
   term1 = 2.0*M*r + Sigma;
   term2 = A - a*a*term1*sin(theta)*sin(theta);

   local_grid->lapse = 1.0/sqrt(1.0 + 2.0*M*r/Sigma);

   local_grid->beta_con[0] = (2.0*M*r/Sigma)/(1.0 + 2*M*r/Sigma);
   local_grid->beta_con[1] = 0.0;
   local_grid->beta_con[2] = 0.0;

   local_grid->gamma_con[0][0] = A*Sigma/(term2*term1);
   local_grid->gamma_con[0][1] = 0.0;
   local_grid->gamma_con[0][2] = Sigma*a/(term2);
   local_grid->gamma_con[1][0] = 0.0;
   local_grid->gamma_con[1][1] = 1.0/Sigma;
   local_grid->gamma_con[1][2] = 0.0;
   local_grid->gamma_con[2][0] = Sigma*a/(term2);
   local_grid->gamma_con[2][1] = 0.0;
   local_grid->gamma_con[2][2] = Sigma/(term2*sin(theta)*sin(theta));

   local_grid->dety = sqrt(term2*term1*sin(theta)*sin(theta)/Sigma);

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
   double M = Black_Hole_Mass;

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
   double phi   = local_grid->x[3];
   double Sigma, Delta, a, M, A;
   double term1, term2, term3;

   a     = Black_Hole_Spin;
   M     = Black_Hole_Mass;
   Sigma = r*r + a*a*cos(theta)*cos(theta);
   Delta = r*r - 2.0*M*r + a*a;
   A     = Sigma*Sigma + a*a*sin(theta)*sin(theta)*(Sigma + 2.0*M*r);
   term1 = 2.0*M*r + Sigma;
   term2 = A - a*a*term1*sin(theta)*sin(theta);
   term3 = 1 + 2.0*M*r/Sigma;

   der->dlapse[0] = (r*r - a*a*cos(theta)*cos(theta))*M/(sqrt(Sigma)*pow(term1,3.0/2.0));
   der->dlapse[1] = -M*a*a*r*sin(2.0*theta)/(Sigma*Sigma*pow(term3,3.0/2.0));
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

   der->dgam[0][0][0] = 2.0*M*(Sigma - 2.0*r*r)/(Sigma*Sigma);
   der->dgam[0][0][1] = 0.0;
   der->dgam[0][0][2] = 2.0*a*M*(2.0*r*r - Sigma)*sin(theta)*sin(theta)/(Sigma*Sigma);
   der->dgam[0][1][0] = 0.0;
   der->dgam[0][1][1] = 2.0*r;
   der->dgam[0][1][2] = 0.0;
   der->dgam[0][2][0] = 2.0*a*M*(2.0*r*r - Sigma)*sin(theta)*sin(theta)/(Sigma*Sigma);
   der->dgam[0][2][1] = 0.0;
   der->dgam[0][2][2] = 2.0*M*a*a*pow(sin(theta),4.0)/Sigma - 4.0*M*a*a*r*r*pow(sin(theta),4.0)/(Sigma*Sigma) + 2.0*r*sin(theta)*sin(theta);

   der->dgam[1][0][0] = 2.0*M*a*a*r*sin(2.0*theta)/(Sigma*Sigma);
   der->dgam[1][0][1] = 0.0;
   der->dgam[1][0][2] = -a*(2.0*M*Sigma*r + 2.0*M*a*a*r*sin(theta)*sin(theta) + Sigma*Sigma)*sin(2.0*theta);
   der->dgam[1][1][0] = 0.0;
   der->dgam[1][1][1] = -a*a*sin(2.0*theta);
   der->dgam[1][1][2] = 0.0;
   der->dgam[1][2][0] = -a*(2.0*M*Sigma*r + 2.0*M*a*a*r*sin(theta)*sin(theta) + Sigma*Sigma)*sin(2.0*theta);
   der->dgam[1][2][1] = 0.0;
   der->dgam[1][2][2] = (4.0*M*Sigma*a*a*r*sin(theta)*sin(theta) + 2.0*M*a*a*a*a*r*pow(sin(theta),4.0) + pow(Sigma,3.0) + Sigma*Sigma*a*a*sin(theta)*sin(theta))*sin(2.0*theta)/(Sigma*Sigma);

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
