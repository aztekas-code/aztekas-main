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
            if(strcmp(t_key,"density_0")==0)
            {
               density_0 = atof(t_value);
            }

            if(strcmp(t_key,"pressure_0")==0)
            {
               pressure_0 = atof(t_value);
            }

            if(strcmp(t_key,"Rad")==0)
            {
               Rad = atof(t_value);
            }

            if(strcmp(t_key,"r_acc")==0)
            {
               r_acc = atof(t_value);
            }

            if(strcmp(t_key,"K")==0)
            {
                if(strcmp(t_value,"1")==0)
                {
                    K = 1.0000000001;
                }
                else if(strcmp(t_value,"2")==0)
                {
                    K = 4./3.;
                }
                else if(strcmp(t_value,"3")==0)
                {
                    K = 5./3.;
                }
                else 
                {
                    K = atof(t_value);
                }                               
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
   
   //allocate space for gas cloud mass vector
   mass = (double *)malloc((Nx1+1)*sizeof(double));

   return 0;   
}

