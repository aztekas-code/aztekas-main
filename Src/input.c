/**
 * @file input.c
 *
 * @author Emilio Tejeda
 *
 * @brief Important input parameters for \a aztekas.
 */

//Do not erase any of these libraries//
#include"main.h"

FILE *paramfile;

int Read_Parameters_File(char const *paramfile_name)
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
            if(strcmp(t_key,"outputdirectory")==0)
            {
               strcpy(outputdirectory,t_value);
            }
         
            if(strcmp(t_key,"outputfile")==0)
            {
               strcpy(outputfile,t_value);
            }

            if(strcmp(t_key,"restart_simulation")==0)
            {
               restart_simulation = atoi(t_value);
            }

            if(strcmp(t_key,"restartfile")==0)
            {
               strcpy(restartfile,t_value);
            }

            if(strcmp(t_key,"restart_filecount")==0)
            {
               restart_filecount = atoi(t_value);
            }

            if(strcmp(t_key,"binary")==0)
            {
               binary = atoi(t_value);
            }

            if(strcmp(t_key,"tmax")==0)
            {
               tmax = atof(t_value);
            }

            if(strcmp(t_key,"timefile")==0)
            {
               timefile = atof(t_value);
            }

            if(strcmp(t_key,"cou")==0)
            {
               cou = atof(t_value);
            }

            if(strcmp(t_key,"K")==0)
            {
               K = atof(t_value);
            }

            if(strcmp(t_key,"Nx1")==0)
            {
               Nx1 = atoi(t_value);
            }

            if(strcmp(t_key,"Nx2")==0)
            {
               Nx2 = atoi(t_value);
            }

            if(strcmp(t_key,"Nx3")==0)
            {
               Nx3 = atoi(t_value);
            }

            if(strcmp(t_key,"x1max")==0)
            {
               x1max = atof(t_value);
            }

            if(strcmp(t_key,"x2max")==0)
            {
               #if alfa == CARTESIAN || alfa == CYLINDRICAL
               x2max = atof(t_value);
               #elif alfa == SPHERICAL
               x2max = M_PI*atof(t_value);
               #endif
            }

            if(strcmp(t_key,"x3max")==0)
            {
               #if alfa == CARTESIAN 
               x3max = atof(t_value);
               #elif alfa == CYLINDRICAL || alfa == SPHERICAL
               x3max = M_PI*atof(t_value);
               #endif
            }

            if(strcmp(t_key,"x1min")==0)
            {
               x1min = atof(t_value);
            }

            if(strcmp(t_key,"x2min")==0)
            {
               #if alfa == CARTESIAN || alfa == CYLINDRICAL
               x2min = atof(t_value);
               #elif alfa == SPHERICAL
               x2min = M_PI*atof(t_value);
               #endif
            }

            if(strcmp(t_key,"x3min")==0)
            {
               x3min = atof(t_value);
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

