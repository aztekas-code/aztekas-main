/**
 * @file /common/print_values.c
 *
 * @authors Alejandro Aguayo-Oritz
 *
 * @brief Print values.
 */

#include"main.h"

void Print_Values()
{
   if(graf == 1)
   {
      if(binary == TRUE)
      {
         Output1_bin(TRUE);
      }
      else
      {
         Output1(TRUE);
      }
   }
   else if(graf == 2)
   {
      if(binary == TRUE)
      {
         Output2_bin(TRUE);
      }
      else
      {
         Output2(TRUE);
      }
   }
   else if(graf == 3)
   {
      if(binary == TRUE)
      {
         Output3_bin(TRUE);
      }
      else
      {
         Output3(TRUE);
      }
   }

   /**
    * Forcing termination
    */
   Alternative_Termination();
}
