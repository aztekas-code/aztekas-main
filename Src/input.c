/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"./Headers/main.h"

FILE *paramfile;

int read_parameters_file(char const *paramfile_name)
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
#if alfa == 0 || alfa == 1
               x2max = atof(t_value);
#elif alfa == 2
               x2max = M_PI*atof(t_value);
#endif
            }

            if(strcmp(t_key,"x3max")==0)
            {
#if alfa == 0 
               x3max = atof(t_value);
#elif alfa == 1 || alfa == 2
               x3max = M_PI*atof(t_value);
#endif
            }

            if(strcmp(t_key,"x1min")==0)
            {
               x1min = atof(t_value);
            }

            if(strcmp(t_key,"x2min")==0)
            {
#if alfa == 0 || alfa == 1
               x2min = atof(t_value);
#elif alfa == 2
               x2min = M_PI*atof(t_value);
#endif
            }

            if(strcmp(t_key,"x3min")==0)
            {
               x3min = atof(t_value);
            }

            if(strcmp(t_key,"r_out")==0)
            {
               r_out = atof(t_value);
            }

            if(strcmp(t_key,"r_in")==0)
            {
               r_in = atof(t_value);
            }

            if(strcmp(t_key,"rho_boundary")==0)
            {
               rho_boundary = atoi(t_value);
            }

            if(strcmp(t_key,"theta_0")==0)
            {
               theta_0 = atof(t_value)*M_PI;
            }

            if(strcmp(t_key,"delta_theta")==0)
            {
               delta_theta = atof(t_value);
            }

            if(strcmp(t_key,"density_0")==0)
            {
               density_0 = atof(t_value);
            }

            if(strcmp(t_key,"pressure_0")==0)
            {
               pressure_0 = atof(t_value);
            }

            if(strcmp(t_key,"velocity_0")==0)
            {
               velocity_0 = atof(t_value);
            }

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

            if(strcmp(t_key,"n_jet")==0)
            {
               n_jet = atof(t_value);
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

            if(strcmp(t_key,"n_atm")==0)
            {
               n_atm = atof(t_value);
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

