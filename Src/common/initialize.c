/*!
 * @file initialize.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Initialize \f$ U \f$ from Restart_File or from initial values
 * defined at the User Directory file: initial.c
 */

#include"main.h"

void Init_Simulation(double *tprint, int *itprint)
{
   if(restart_simulation == TRUE)
   {
      if(access(restartfile,F_OK) != -1)
      {
         printf("Restarting from file ");
         for (int c = 0; restartfile[c] != '\0'; c++)
         {
            printf("%c",restartfile[c]);
         }
         printf("\n");
         printf("\n");
      }
      else
      {
         printf("Restarting file ");
         for (int c = 0; restartfile[c] != '\0'; c++)
         {
            printf("%c",restartfile[c]);
         }
         printf(" does not exist.\n");
         printf("\n");
         exit(EXIT_FAILURE);
      }


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
