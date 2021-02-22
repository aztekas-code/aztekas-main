/**
 * @file user_input.c
 *
 * @author Emilio Tejeda
 *
 * @brief Important input parameters for \a aztekas.
 */

//Do not erase any of these libraries//
#include"main.h"

FILE *paramfile;

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

            if(strcmp(t_key,"Black_Hole_Spin")==0)
            {
               Black_Hole_Spin = atof(t_value);
            }

            if(strcmp(t_key,"cs")==0)
            {
               cs = atof(t_value);
            }

            if(strcmp(t_key,"Temp")==0)
            {
               Temp = atof(t_value);
            }

            if(strcmp(t_key,"r_star")==0)
            {
               r_star = atof(t_value);
            }

            if(strcmp(t_key,"density_0")==0)
            {
               density_0 = atof(t_value);
            }

            if(strcmp(t_key,"velocity_0")==0)
            {
               velocity_0= atof(t_value);
            }

            if(strcmp(t_key,"last")==0)
            {
               strcpy(last,t_value);
            }

            if(strcmp(t_key,"MDOT_ERR")==0)
            {
               MDOT_ERR = atof(t_value);
            }

            if(strcmp(t_key,"EoS")==0)
            {
               strcpy(eqstate,t_value);
            }
         }
      }
            
      // Ignore the rest of the line.
      int fscanret = fscanf(paramfile, "%*[^\n]"); 
   }
      
   fclose(paramfile);         
   
   /**
    * Get h and c_s from Temperature
    */
   eos_ eos;
   gauge_ local_grid;
   double P[2];
   P[0] = 1.0;
   P[1] = Temp;

   EoS(&eos,P,&local_grid);

   /**
    * Get equivalent polytropic index K_equiv
    */
   double K_h = (eos.h - 1.0)/Temp; // For ideal-polytropic gas, K_h = K/(K - 1).
   double K_equiv = K_h/(K_h - 1.0);// For ideal-polytropic gas, K_equiv = K.
   K = K_equiv;
   
   /**
    * Use c_s to parametrize the size of the domain x1max
    */
   cs = eos.cs;
   if(strcmp(eqstate,"Stiff") == 0)
   {
      K = 1.333333333;
      cs -= 0.001;
   }

   rB = 1.0/pow(cs,2.0);

   if(40*rs(K,cs) >= 10*rB)
   {
      x1max = 40*rs(K,cs);
      x1min = x1min*rs(K,cs);
   }
   else
   {
      x1max = x1max*rB;
      x1min = x1min*rs(K,cs);
   }

   x1max = 300.0;

   if(x1min < 1.5)
   {
      x1min = 1.5;
   }

   /**
    * Compute max time in terms of c_s
    */
   tB       = x1max/cs;
   tmax     = tmax*tB;
   timefile = timefile*tB;

   printf("Sonic radius r_s = %f\n",rs(K,cs));
   printf("Bondi radius r_B = %f\n",rB);
   printf("Temperature = %.4e\n",Temp);
   printf("Speed of sound = %f\n",cs);
   printf("Angular momentum = %f\n",l_0);
   printf("Equivalent polytropic index = %f\n",K);
   printf("Rmax = %f M, Rmax = %f rs, Rmax = %f rB\n",x1max,x1max/rs(K,cs),x1max/rB);
   printf("Rmin = %f M, Rmin = %f rs, Rmin = %f rB\n",x1min,x1min/rs(K,cs),x1min/rB);
   printf("Max time = %e M, %f tB\n",tmax,tmax/tB);

   if(rs(K_equiv,cs) > x1max)
   {
      printf("######################################\n");
      printf("# Sonic radius is larger than domain #\n");
      printf("######################################\n");
      exit(EXIT_FAILURE);
   }

   if(*outputfile == '\0')
   {
      printf("Check parameters file: Didn't find a name for output file\n");
      exit(EXIT_FAILURE);
   }

   MDOT_TIME = pow(10,floor(log10(tmax)) - floor(log10(MDOT_DATA)));

   return 0;   
}
