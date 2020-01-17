/*
 * File Name : user_input.c
 * Description : aztekas user input parameters for Relativistic Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:41:05
 * Created By : Alejandro Aguayo-Ortiz
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
            if(strcmp(t_key,"rhotl")==0)
            {
               rhotl = atof(t_value);
            }

            if(strcmp(t_key,"ptl")==0)
            {
               ptl = atof(t_value);
            }

            if(strcmp(t_key,"vx1tl")==0)
            {
               vx1tl = atof(t_value);
            }

            if(strcmp(t_key,"vx2tl")==0)
            {
               vx2tl = atof(t_value);
            }

            if(strcmp(t_key,"rhotr")==0)
            {
               rhotr = atof(t_value);
            }

            if(strcmp(t_key,"ptr")==0)
            {
               ptr = atof(t_value);
            }

            if(strcmp(t_key,"vx1tr")==0)
            {
               vx1tr = atof(t_value);
            }

            if(strcmp(t_key,"vx2tr")==0)
            {
               vx2tr = atof(t_value);
            }

            if(strcmp(t_key,"rhobl")==0)
            {
               rhobl = atof(t_value);
            }

            if(strcmp(t_key,"pbl")==0)
            {
               pbl = atof(t_value);
            }

            if(strcmp(t_key,"vx1bl")==0)
            {
               vx1bl = atof(t_value);
            }

            if(strcmp(t_key,"vx2bl")==0)
            {
               vx2bl = atof(t_value);
            }

            if(strcmp(t_key,"rhobr")==0)
            {
               rhobr = atof(t_value);
            }

            if(strcmp(t_key,"pbr")==0)
            {
               pbr = atof(t_value);
            }

            if(strcmp(t_key,"vx1br")==0)
            {
               vx1br = atof(t_value);
            }

            if(strcmp(t_key,"vx2br")==0)
            {
               vx2br = atof(t_value);
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

