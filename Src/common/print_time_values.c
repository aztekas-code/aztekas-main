/**
 * @file /common/print_time_values.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Print values each time.
 */

#include"main.h"

void Print_Time_Values(double *tprint, double *dtprint, int it)
{
   /**
    * Print the Mass Accretion Rate \f[ \dot{M} \f]
    */
   #if MDOT == TRUE
   Mass_Accretion_Rate(U);
   #endif

   if(grid.time >= *tprint || CHECK_NAN == TRUE)
   {
      printf("Time = %e, dt = %e\n",grid.time,dt);
//      if(graf == 1)
//      {
         if(binary == TRUE)
         {
            Output_bin(it);
         }
         else
         {
            Output_ascii(it);
         }
//      }
//      else if(graf == 2)
//      {
//         if(binary == TRUE)
//         {
//            Output2_bin(it);
//         }
//         else
//         {
//            Output2(it);
//         }
//      }
//      else if(graf == 3)
//      {
//         if(binary == TRUE)
//         {
//            Output3_bin(it);
//         }
//         else
//         {
//            Output3(it);
//         }
//      }

      /**
       * Increase the time of printing by dtprint and increase
       * the humber of the output file
       */
      *tprint = *tprint + *dtprint;
      itprint = it + 1;

      /**
       * Forcing termination
       */
      Alternative_Termination();
   }
}
