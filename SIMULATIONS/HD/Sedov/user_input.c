/*
 * File Name : user_input.c
 * Description : aztekas user input parameters for Sedov Blast Wave
 * Creation Date : 26-09-2019
 * Last Modified : 26-09-2019 23:58:14
 * Created By : Emilio Tejeda
 */

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
            if(strcmp(t_key,"n_0")==0)
            {
               rho_0 = atof(t_value);
            }

            if(strcmp(t_key,"p_0")==0)
            {
               p_0 = atof(t_value);
            }

            if(strcmp(t_key,"E_0")==0)
            {
               E_0 = atof(t_value);
            }

            if(strcmp(t_key,"x_0")==0)
            {
               x_0 = atof(t_value);
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
