/**
 * @mainpage AZTEKAS: a hydrodynamic GPL code
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file main.c
 *
 * @author Alejandro Aguayo-Ortiz.
 *
 * @brief Main file of aztekas.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<omp.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"
#include"param.h"

int main(int argc, char* argv[])
{
   int itprint;
   double dtprint, tprint;

   if(argc != 2) 
   {
      printf("%s\n","Wrong number of arguments") ;
      printf("%s\n","Execute as:") ;    
      printf("%s\n","./aztekas paramfile") ;  
      exit(EXIT_FAILURE);
   }

   strcpy(paramfile_name, argv[1]);
   
   read_parameters_file(paramfile_name);   
      
	// create output directory
   char create_dir[] = "mkdir -p ";	
	strcat(create_dir,outputdirectory);	
	int sysret = system(create_dir);

   // Include the ghost cells
   new_SIZE();

   allocateArray();

   //We set the mesh func_planarMESH.c
   MESH();

   //Time interval between data dumps
   dtprint = timefile;

   //We set the initial parameters func_planarINITIAL.c
   if( restart_simulation == 1 )
   {
      RESTART();

      tprint = time;
      itprint = restart_filecount;
   }
   else
   {
      INITIAL();
      tprint  = 0.0; //Initialize printing parameter
      itprint = 0;   //Initialize file numeration
   }

   start = omp_get_wtime();
   while(time <= tmax)
   {
      //In this part we compute the time step
      dt = TIMESTEP();

      //We print the values: file (DATOS*) and to terminal func_planarOUTPUT.c
      PrintValues(&tprint,&dtprint,&itprint);

      //In here we set the integration method (Finite volume method)
      INTEGRATION();
   }

   PrintValues(&tprint,&dtprint,&itprint);

   delta = omp_get_wtime() - start;
   printf("Delta %.4g seconds with %d threads\n",delta,4);

   free(X1);
   free(X2);
   free(X3);

   return 0;
}
