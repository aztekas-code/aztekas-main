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
            if(strcmp(t_key,"rho_jet")==0)
            {
               rho_jet = atof(t_value);
            }

            if(strcmp(t_key,"p_jet")==0)
            {
               p_jet = atof(t_value);
            }

            if(strcmp(t_key,"vx1_jet")==0)
            {
               vx1_jet = atof(t_value);
            }

            if(strcmp(t_key,"vx2_jet")==0)
            {
               vx2_jet = atof(t_value);
            }

            if(strcmp(t_key,"vx3_jet")==0)
            {
               vx3_jet = atof(t_value);
            }

            if(strcmp(t_key,"rho_atm")==0)
            {
               rho_atm = atof(t_value);
            }

            if(strcmp(t_key,"p_atm")==0)
            {
               p_atm = atof(t_value);
            }

            if(strcmp(t_key,"vx1_atm")==0)
            {
               vx1_atm = atof(t_value);
            }

            if(strcmp(t_key,"vx2_atm")==0)
            {
               vx2_atm = atof(t_value);
            }

            if(strcmp(t_key,"vx3_atm")==0)
            {
               vx3_atm = atof(t_value);
            }

            if(strcmp(t_key,"r_jet")==0)
            {
               r_jet = atof(t_value);
            }

            if(strcmp(t_key,"z_jet")==0)
            {
               z_jet = atof(t_value);
            }
         }
      }
            
      // Ignore the rest of the line.
      int fscanret = fscanf(paramfile, "%*[^\n]"); 
   }
      
   fclose(paramfile);         

   if(*outputfile == '\0')
   {
      printf("Check parameters file: Didn't find a name for output file\n");
      exit(EXIT_FAILURE);
   }

   return 0;   
}

