#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define lfac   1.0
#define MAX(a,b) (((a)<(b))?(a):(b))
FILE *paramfile;
int N;
double c_infty, Gamma, Temp, log_Temp, dT, K;
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
   double M    = 1.0;

   file = fopen("mdot.dat","w");

   log_Temp  = -5;
   Gamma = 1.001;
   dT    = (3.0 + 5.0)/1000.0;
   rho_infty = 1.0;
   strcpy(EoS,"Ideal");
   while(log_Temp <= 3)
   {
      Temp = pow(10.0,log_Temp);

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

      double Mdot = 4*M_PI*rho_star*U_star*r_star*r_star/(M_PI*rho_infty);

      fprintf(file,"%e %e\n",Temp,-Mdot);

      log_Temp = log_Temp + dT;
   }

   fclose(file);


   return 0;
}
