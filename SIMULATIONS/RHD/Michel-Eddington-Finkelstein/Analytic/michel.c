#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define MAX(a,b) (((a)<(b))?(a):(b))

int main(int argc, char *argv[])
{
   int N = 128;
   int gc = 3;
   // Define polytropic index Gamma, rho_infty, c_infty and M.
   // xmax and xmin
   double M         = 1.0;
   double Gamma     = 5.0/3.0;
   double rho_infty = 1.0;
   double c_infty   = 0.1;
   double xmin = 1.5;
   double xmax = 40;

   // Compute the polytropic constant K
   double K;

   K = (Gamma - 1.0)*pow(c_infty,2.0)/((Gamma - 1.0 - pow(c_infty,2.0))*(Gamma*pow(rho_infty,Gamma - 1.0)));

   // Compute speed of sound at the sonic point c_star
   double c_star;
   double c_star2;
   double arg1, arg2;
   double numerator, denominator;
   
   numerator   = 54.0*pow(Gamma,3.0) - 351.0*pow(Gamma,2.0) + 558.0*Gamma + \
                 486.0*(Gamma - 1.0)*pow(c_infty,2.0) - 243.0*pow(c_infty,4.0) - \
                 259.0;
   denominator = 2.0*pow(3.0*Gamma - 2.0,3.0);
   arg1 = numerator/denominator;
   arg2 = M_PI/3.0 + acos(arg1)/3.0;
   c_star2 = (6.0*Gamma - 7.0 + 2*(3.0*Gamma - 2.0)*cos(arg2))/9.0;
   c_star  = fabs(sqrt(c_star2));

   // Compute the sonic point r_star
   double r_star;

   r_star = M*(1.0 + 3.0*c_star2)/(2.0*c_star2);

   // Compute the radial component of the four velocity at the sonic point
   double U_star;

   U_star = - sqrt(c_star2/(1.0 + 3.0*c_star2));

   // Compute energy density at the sonic point
   double e_star;

   e_star = c_star2/(Gamma*(Gamma - 1.0 - c_star2));

   // Compute enthalpy at the sonic point
   double h_star;

   h_star = 1.0 + Gamma*e_star;

   // Compute the density at the sonic point
   double rho_star;

   rho_star = pow((Gamma - 1.0)*e_star/K,1.0/(Gamma - 1.0));

   // Compute constant B = rho_star*r_star²*U_star
   double B = r_star*r_star*rho_star*U_star;

   //Print values
   printf("############################################\n");
   printf("############## Michel Problem ##############\n");
   printf("############################################\n");
   printf("# \n");
   printf("# This is the solution of the Michel Problem\n");
   printf("# with parameters\n");
   printf("# \n");
   printf("# Gamma     = %f\n",Gamma);
   printf("# c_infty   = %e\n",c_infty);
   printf("# rho_infty = %e\n",rho_infty);
   printf("# M = %e\n",M);
   printf("# r_min = %e and r_max = %e\n",xmin,xmax);
   printf("# \n");
   printf("# from which we obtain\n");
   printf("# \n");
   printf("# c_star   = %e\n",c_star);
   printf("# r_star   = %e\n",r_star);
   printf("# rho_star = %e\n",r_star);
   printf("# U_star   = %e\n",U_star);
   printf("# B        = %e\n",B);
   printf("# K        = %e\n",K);
   printf("# \n");
   printf("############################################\n");

   // Define de grid mesh
   double X[N+1+2*gc];
   double dx = (xmax - xmin)/((double)N);
   if( strcmp(argv[1],"log") == 0)
   {
      for(int i = 0; i <= N+2*gc; i++)
      {
         X[i] = xmin + exp(log((xmax - xmin + 1.0))*(i-gc)/(N)) - 1.0;
      }
   }
   else
   {
      for(int i = 0; i <= N+2*gc; i++)
      {
         X[i] = xmin + (i-gc)*dx;
      }
   }

   // Compute the radial component of the four velocity at each point 
   int count, count_2=0;
   int i, j;
   double RHO[N+2*gc+1], PRE[N+2*gc+1], U[N+2*gc+1], V[N+2*gc+1];
   double r, rho, a;
   double f, df, ddf, u0, u, u_last;
   double E, C;
   double p;

   u0 = sqrt(1.0 - 2.0*M/xmin + 3.0);
   i  = 0;

   while (r < r_star)
   {
      f = 1.0;
      r = X[i];

      E = ((Gamma - 1.0)/(K*Gamma))*((pow(r*r,Gamma - 1.0))/pow(fabs(B),Gamma - 1.0)); 
      C = 1.0 - pow(c_infty,2.0)/(Gamma - 1.0);

      a = (Gamma - 1.0);

      LOOP:do{
         
         f  = 1.0 - 1.0/(1.0 + E*pow(u0,a)) - sqrt(1.0 - 2.0*M/r + u0*u0)*C;
         df = E*a*pow(u0,a-1.0)/pow(1.0 + E*pow(u0,a),2.0) - \
              u0*C/(sqrt(1.0 - 2.0*M/r + u0*u0));

         u = u0 - f/df;
         u_last = u0;
         u0 = u;

         rho = -B/(r*r*u);

      }while (fabs(f) > 0.00000001);

      RHO[i] = rho;
      PRE[i] = K*pow(rho,Gamma);
      U[i]   = -u;
      double lapse = sqrt(1.0/(1.0 + 2.0*M/r));
      double br    = (2.0*M/r)*(1.0/(1.0 + 2.0*M/r));
      double gtt   = - (1.0 - 2.0*M/r);
      double gTT   = - 1.0 - 2.0*M/r;
      double gtr   = 2.0*M/r;
      double gTR   = 2.0*M/r;
      double grr   = 1.0 + 2.0*M/r;
      double gRR   = 1.0 - 2.0*M/r;

      double Uur   = U[i];
      double Udt   = - sqrt(Uur*Uur - gtt);
      double Udr   = (gtr*gTT*Udt + grr*Uur)/(1.0 - gtr*gTR);
      double Uut   = gTT*Udt + gTR*Udr;

      double W     = lapse*Uut;

      V[i] = (Udr / W);

      j = i;
      i++;
   }

   u0 = 0.001;

   for(int i = N+2*gc; i >= j; i--)
   {
      f = 1.0;
      r = X[i];

      E = ((Gamma - 1.0)/(K*Gamma))*((pow(r*r,Gamma - 1.0))/pow(fabs(B),Gamma - 1.0)); 
      C = 1.0 - pow(c_infty,2.0)/(Gamma - 1.0);

      a = (Gamma - 1.0);

      do{
         f  = 1.0 - 1.0/(1.0 + E*pow(u0,a)) - sqrt(1.0 - 2.0*M/r + u0*u0)*C;
         df = E*a*pow(u0,a-1.0)/pow(1.0 + E*pow(u0,a),2.0) - \
              u0*C/(sqrt(1.0 - 2.0*M/r + u0*u0));

         u = u0 - f/df;
         u_last = u0;
         u0 = u;

         rho = -B/(r*r*u);

      }while (fabs(f) > 0.00000001);

      RHO[i] = rho;
      PRE[i] = K*pow(rho,Gamma);
      U[i]   = -u;

      double lapse = sqrt(1.0/(1.0 + 2.0*M/r));
      double br    = (2.0*M/r)*(1.0/(1.0 + 2.0*M/r));
      double gtt   = - (1.0 - 2.0*M/r);
      double gTT   = - 1.0 - 2.0*M/r;
      double gtr   = 2.0*M/r;
      double gTR   = 2.0*M/r;
      double grr   = 1.0 + 2.0*M/r;
      double gRR   = 1.0 - 2.0*M/r;

      double Uur   = U[i];
      double Udt   = - sqrt(Uur*Uur - gtt);
      double Udr   = (gtr*gTT*Udt + grr*Uur)/(1.0 - gtr*gTR);
      double Uut   = gTT*Udt + gTR*Udr;

      double W     = lapse*Uut;

      V[i] = (Udr / W);
   }

   for(int i = 0; i <= N+2*gc; i++)
      printf("%e %e %e %e\n",X[i],RHO[i],PRE[i],V[i]);
}
