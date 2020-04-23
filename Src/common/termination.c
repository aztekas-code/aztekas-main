/**
 * @file termination.c
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Termination functions
 */

#include"main.h"

void Termination()
{
   Ending_Message();
}

void Alternative_Termination()
{
#if MDOT == TRUE                                                          
   if(Mdot_end == TRUE)                                                         
   {   
      printf("\n");
      printf("The mean value of Mdot remain constant in %e\n",MDOT_ERR);
      Ending_Message();
      exit(EXIT_FAILURE);
   }
#endif
}
