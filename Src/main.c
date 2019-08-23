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
#include"main.h"

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
   
   Read_Parameters_File(paramfile_name);   
   User_Parameters(paramfile_name);   
      
	// create output directory
   char create_dir[] = "mkdir -p ";	
	strcat(create_dir,outputdirectory);	
	int sysret = system(create_dir);

   // Print Simulation Parameters
   CheckSimParameters();
   if (check_param == TRUE) getchar();

   // Include the ghost cells
   New_Size();

   Allocate_Array();

   //We set the mesh func_planarMESH.c
   Mesh();

   //Time interval between data dumps
   dtprint = timefile;

   //We set the initial parameters func_planarINITIAL.c
   if( restart_simulation == TRUE )
   {
      if( binary == TRUE )
      {
         Restart_Bin();
      }
      else
      {
         Restart();
      }

      tprint = grid.time;
      itprint = restart_filecount;
   }
   else
   {
      Initial();
      tprint  = 0.0; //Initialize printing parameter
      itprint = 0;   //Initialize file numeration
      U0 = U;
   }

   start = omp_get_wtime();
   while(grid.time <= tmax)
   {
      //In this part we compute the time step
      dt = TimeStep();

      //We print the values: file (DATOS*) and to terminal func_planarOUTPUT.c
      PrintValues(&tprint,&dtprint,&itprint);

      //In here we set the integration method (Finite volume method)
      Integration();

      printf("Time = %e, dt = %e\r",grid.time,dt);
      fflush(stdout); 
   }

   PrintValues(&tprint,&dtprint,&itprint);

   delta = omp_get_wtime() - start;
   printf("Delta %.4g seconds with %d threads\n",delta,1);

   free(grid.X1);
   free(grid.X2);
   free(grid.X3);

   return 0;
}
