/*
 * File Name : user_input.c
 * Description : aztekas user input parameters for Rayleigh Taylor
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:02:06
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

      if(t_firstChar != '%')
      {
         // Not a comment so read the key value pair
         // Move back one space in the input stream with seek
         fseek(paramfile, -1, SEEK_CUR);

         if(fscanf (paramfile, "%s = %s", t_key, t_value) == 2)
         {
            if(strcmp(t_key,"rhod")==0)
            {
               rhod = atof(t_value);
            }

            if(strcmp(t_key,"pd")==0)
            {
               pd = atof(t_value);
            }

            if(strcmp(t_key,"vx1d")==0)
            {
               vx1d = atof(t_value);
            }

            if(strcmp(t_key,"vx2d")==0)
            {
               vx2d = atof(t_value);
            }

            if(strcmp(t_key,"vx3d")==0)
            {
               vx3d = atof(t_value);
            }

            if(strcmp(t_key,"rhou")==0)
            {
               rhou = atof(t_value);
            }

            if(strcmp(t_key,"pu")==0)
            {
               pu = atof(t_value);
            }

            if(strcmp(t_key,"vx1u")==0)
            {
               vx1u = atof(t_value);
            }

            if(strcmp(t_key,"vx2u")==0)
            {
               vx2u = atof(t_value);
            }

            if(strcmp(t_key,"vx3u")==0)
            {
               vx3u = atof(t_value);
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

