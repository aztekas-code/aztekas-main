#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define lfac   1.0
#define MAX(a,b) (((a)<(b))?(a):(b))
FILE *paramfile;
int N;
double c_infty, Gamma, Temp;
double rho_infty;
double x1min, x1max;
char EoS[50];

double rs(double K, double c_s)
{
   double n       = 1.0/(K - 1.0);
   double h_inf   = 1.0/(1.0 - n*pow(c_s,2.0));
   double Psi     = acos((3.0/(2.0*n*h_inf))*pow((n + 3.0)/(3.0*n),-1.5))/3.0;
   double h_c     = 2.0*h_inf*sqrt((n + 3.0)/(3.0*n))*sin(Psi + M_PI/6.0);
   double c_s_c2  = (h_c - 1.0)/(n*h_c);
   double v_c2    = c_s_c2/(1.0 + 3.0*c_s_c2);
   double rs      = 0.5/v_c2;

   return rs;
}

int User_Parameters(char const *paramfile_name)
{
   int   BUFFER_SIZE = 612;
   char t_key[BUFFER_SIZE], t_value[BUFFER_SIZE], t_firstChar;
   
   //if input paramfile_name is null
   if(!paramfile_name || *paramfile_name == '\0')
   {
      fprintf(stderr, "Invalid parameter file name\n");
      exit(EXIT_FAILURE);
   }
      
   //open file
   paramfile = fopen (paramfile_name, "r");

   if (paramfile == NULL)
   {
      fprintf(stderr, "Error opening parameter file: %s\n", paramfile_name);
      exit(EXIT_FAILURE);
   }
      
   //read individual settings lines
   while(fscanf(paramfile, " %c", &t_firstChar) == 1 )
   {

      if(t_firstChar != '/' && t_firstChar != '%')
      {
         // Not a comment so read the key value pair
         // Move back one space in the input stream with seek
         fseek(paramfile, -1, SEEK_CUR);

         if(fscanf (paramfile, "%s = %s", t_key, t_value) == 2)
         {

            if(strcmp(t_key,"Nx1")==0)
            {
               N = atoi(t_value);
            }

            if(strcmp(t_key,"cs")==0)
            {
               c_infty = atof(t_value);
            }

            if(strcmp(t_key,"Temp")==0)
            {
               Temp = atof(t_value);
            }

            if(strcmp(t_key,"K")==0)
            {
               Gamma = atof(t_value);
            }

            if(strcmp(t_key,"x1max")==0)
            {
               x1max = atof(t_value);
            }

            if(strcmp(t_key,"x1min")==0)
            {
               x1min = atof(t_value);
            }

            if(strcmp(t_key,"EoS")==0)
            {
               strcpy(EoS,t_value);
            }
         }
      }
            
      // Ignore the rest of the line.
      int fscanret = fscanf(paramfile, "%*[^\n]"); 
   }

   if(strcmp(EoS,"Ideal") == 0)
   {
      c_infty = sqrt(Gamma*Temp*(Gamma - 1.0)/(Gamma*Temp + Gamma - 1.0));
   }
   else if(strcmp(EoS,"Ryu") == 0)
   {
      double h = 2.0*(6.0*Temp*Temp + 4.0*Temp + 1.0)/(3.0*Temp + 2.0);
      c_infty = sqrt((Temp*(3.0*Temp + 2.0)*(18.0*Temp*Temp + 24.0*Temp + 5.0))/\
                (3.0*(6.0*Temp*Temp + 4.0*Temp + 1.0)*(9.0*Temp*Temp + 12.0*Temp + 2.0)));

      double K_h = (h - 1.0)/Temp;
      Gamma = K_h/(K_h - 1.0);
   }
   else if (strcmp(EoS,"Stiff") == 0)
   {
      c_infty = sqrt((Gamma - 1.0)) - 0.001;
   }

   double rB = 1.0/pow(c_infty,2.0);

   if(40*rs(Gamma,c_infty) >= 10*rB)
   {
      x1max = 40*rs(Gamma,c_infty);
      x1min = x1min*rs(Gamma,c_infty);
   }
   else
   {
      x1max = x1max*rB;
      x1min = x1min*rs(Gamma,c_infty);
   }

   if(x1min < 1.1)
   {
      x1min = 1.1;
   }
      
   fclose(paramfile);         

   return 0;   
}

int main(int argc, char *argv[])
{
   FILE *file;
   User_Parameters(argv[2]);
   int gc = 3;
   // Define polytropic index Gamma, rho_infty, c_infty and M.
   // xmax and xmin
   double M    = 1.0;
   double xmax = x1max;
   double xmin = x1min;

   // Compute the polytropic constant K
   double K;
   rho_infty = 1.0;

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
//   printf("############################################\n");
//   printf("############## Michel Problem ##############\n");
//   printf("############################################\n");
//   printf("# \n");
//   printf("# This is the solution of the Michel Problem\n");
//   printf("# with parameters\n");
//   printf("# \n");
//   printf("# Gamma     = %f\n",Gamma);
//   printf("# c_infty   = %e\n",c_infty);
//   printf("# rho_infty = %e\n",rho_infty);
//   printf("# M = %e\n",M);
//   printf("# r_min = %e and r_max = %e\n",xmin,xmax);
//   printf("# \n");
//   printf("# from which we obtain\n");
//   printf("# \n");
//   printf("# c_star   = %e\n",c_star);
//   printf("# r_star   = %e\n",r_star);
//   printf("# rho_star = %e\n",rho_star);
//   printf("# U_star   = %e\n",U_star);
//   printf("# B        = %e\n",B);
//   printf("# K        = %e\n",K);
//   printf("# Mdot     = %e\n",4*M_PI*rho_star*U_star*r_star*r_star/(M_PI*rho_infty));
//   printf("# \n");
//   printf("############################################\n");

   // Define de grid mesh
   double X[N+1+2*gc];
   double dx = (xmax - xmin)/((double)N);
   if( strcmp(argv[1],"log") == 0)
   {
      for(int i = 0; i <= N+2*gc; i++)
      {
         X[i] = xmin + lfac*exp(log((xmax - xmin + lfac)/lfac)*(i-gc)/(N)) - lfac;
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

   u0 = sqrt(fabs(1.0 - 2.0/(X[0]-0.1)));
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

         //printf("%e %e %e %e\n",r,E,C,a);
         //printf("%e %e %e\n",u0,f,df);
         //getchar();
         u = u0 - f/df;
         u_last = u0;
         u0 = u;

         rho = -B/(r*r*u);
         //printf("B: %e, r: %e, u: %e, rho: %e\n",B,r,u,rho);
         //getchar();

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

      if(i == N + 2*gc + 1)
         goto END;
   }

   // T > 0.01
   u0 = 0.000000001;

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
      U[i] = Uur + W*br/lapse;
   }
    
   END:

   file = fopen("analytic.dat","w");

   for(int i = 0; i <= N+2*gc; i++)
      fprintf(file,"%e %e %e %e %e\n",X[i],RHO[i]/RHO[N+gc],PRE[i]/RHO[N+gc],V[i],U[i]);
   //for(int i = gc; i <= N+gc; i++)
   //   fprintf(file,"%e %e %e %e\n",X[i],RHO[i],PRE[i],V[i]);

   fclose(file);

   printf("Michel Mass Accretion Rate = %e\n",4*M_PI*rho_star*U_star*r_star*r_star/(rho_infty*M_PI));
   printf("Density %e\n",RHO[N+gc]);
   printf("Velocity %e\n",V[N+gc]);

   return 0;
}
