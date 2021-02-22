/**
 * @file /common/print_values.c
 *
 * @authors Alejandro Aguayo-Oritz
 *
 * @brief Print values.
 */

#include"main.h"

void Print_Values0()
{
   char id[50];
   strcpy(id,"out");

   if(graf == 1)
   {
      if(binary == TRUE)
      {
         Output1_bin(id);
      }
      else
      {
         Output1(id);
      }
   }
   else if(graf == 2)
   {
      if(binary == TRUE)
      {
         Output2_bin(id);
      }
      else
      {
         Output2(id);
      }
   }
   else if(graf == 3)
   {
      if(binary == TRUE)
      {
         Output3_bin(id);
      }
      else
      {
         Output3(id);
      }
   }

   /**
    * Forcing termination
    */
   Alternative_Termination();
}

void Print_Values1(char *file_id)
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

   if(graf == 1)
   {
      if(binary == TRUE)
      {
         Output1_bin(id);
      }
      else
      {
         Output1(id);
      }
   }
   else if(graf == 2)
   {
      if(binary == TRUE)
      {
         Output2_bin(id);
      }
      else
      {
         Output2(id);
      }
   }
   else if(graf == 3)
   {
      if(binary == TRUE)
      {
         Output3_bin(id);
      }
      else
      {
         Output3(id);
      }
   }

   /**
    * Forcing termination
    */
   Alternative_Termination();
}
