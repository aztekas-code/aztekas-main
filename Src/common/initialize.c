/*
 * File Name : initialize.c
 * Description :
 * Creation Date : 22-10-2019
 * Last Modified : 22-10-2019 15:21:50
 * Created By :
 */

#include"main.h"

void Init_Simulation(double *tprint, int *itprint)
{
   if(restart_simulation == TRUE)
   {
      if(binary == TRUE)
      {
         Restart_Bin();
      }
      else
      {
         Restart();
      }

      *tprint  = grid.time;          // Initialize time
      *itprint = restart_filecount;  // Initialize number of files
      
   }
   else
   {
      Initial();
      *tprint  = 0.0;  // Initialize time
      *itprint = 0;    // Initialize number of files
   }

   U0 = U;        // Save U(t_n), initial or previous step vector
   Boundaries(U); // Fill with correct boundary conditions.
}
