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

int main(int argc, char *argv[])
{
   /**
    * Check if AZTEKAS is run properly
    */
   Check_Paramfile(paramfile_name,argc,argv);

   /*
    * Read necessary and user defined parameters from file.param
    * and print info in screen
    */
   Default_Parameters(paramfile_name);
   User_Parameters(paramfile_name);

   /**
    * Check paramfile, print info of the simulation on screen
    * and also inside a directory named INFO inside the simulation
    * directory. Print INFO on the screen.
    */
   Manage_Simulation_Info(argc,argv);

   /**
    * Allocate the space for all the arrays used by aztekas
    */
   Allocate_Array();

   /**
    * Create a Cartesian-like mesh grid
    */
   Mesh();

   /**
    * Initialize solution vector U
    */
   Init_Simulation(&tprint,itprint);

   /*
    * Frequency of printing output
    */
   Frequency_Output(&dtprint);

   /**
    * Starts computational spending time
    */
   Computing_Time_Start();

   /**
    * Solves the set of equations
    */
   Equation_System_Solver();

   /**
    * Print final message with computing time
    */
   Ending_Message();
                                                                                
   free(U);                                                                     
   free(grid.X1);                                                               
   free(grid.X2);                                                               
   free(grid.X3); 

   return 0;
}
