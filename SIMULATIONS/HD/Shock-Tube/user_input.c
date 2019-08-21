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
            if(strcmp(t_key,"nl")==0)
            {
               nl = atof(t_value);
            }

            if(strcmp(t_key,"pl")==0)
            {
               pl = atof(t_value);
            }

            if(strcmp(t_key,"vx1l")==0)
            {
               vx1l = atof(t_value);
            }

            if(strcmp(t_key,"vx2l")==0)
            {
               vx2l = atof(t_value);
            }

            if(strcmp(t_key,"vx3l")==0)
            {
               vx3l = atof(t_value);
            }

            if(strcmp(t_key,"nr")==0)
            {
               nr = atof(t_value);
            }

            if(strcmp(t_key,"pr")==0)
            {
               pr = atof(t_value);
            }

            if(strcmp(t_key,"vx1r")==0)
            {
               vx1r = atof(t_value);
            }

            if(strcmp(t_key,"vx2r")==0)
            {
               vx2r = atof(t_value);
            }

            if(strcmp(t_key,"vx3r")==0)
            {
               vx3r = atof(t_value);
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

