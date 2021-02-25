/**
 * @file /common/print_values.c
 *
 * @authors Alejandro Aguayo-Oritz
 *
 * @brief Print values.
 */

#include"main.h"

void Print_Values_0()
{
   char id[50];
   strcpy(id,"out");

   if(binary == TRUE)
   {
      Output_bin(id);
   }
   else
   {
      Output_ascii(id);
   }

   /**
    * Forcing termination
    */
   Alternative_Termination();
}

void Print_Values_1(char *file_id)
{
   char id[50];
   if (strcmp(file_id,"") == 0)
   {
      strcpy(id,"out");
   }
   else
   {
      strcpy(id,file_id);
   }

   if(binary == TRUE)
   {
      Output_bin(id);
   }
   else
   {
      Output_ascii(id);
   }

   /**
    * Forcing termination
    */
   Alternative_Termination();
}
