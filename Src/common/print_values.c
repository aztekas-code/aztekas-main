/**
 * @file /common/print_values.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Print values.
 */

#include"main.h"

void Print_Values(double *tprint, double *dtprint, int *itprint)
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
      if(graf == 1)
      {
         if(binary == 1)
         {
            Output1_bin(itprint);
         }
         else
         {
            Output1(itprint);
         }
      }
      else if(graf == 2)
      {
         if(binary == 1)
         {
            Output2_bin(itprint);
         }
         else
         {
            Output2(itprint);
         }
      }
      else if(graf == 3)
      {
         if(binary == 1)
         {
            Output3_bin(itprint);
         }
         else
         {
            Output3(itprint);
         }
      }

      /**
       * Increase the time of printing by dtprint and increase
       * the humber of the output file
       */
      *tprint = *tprint + *dtprint;
      ++*itprint;

      #if MDOT == TRUE
      if(Mdot_end == TRUE)
      {                                                                            
         printf("\n");                                                             
         printf("AZTEKAS termination\n");                                          
         printf("The mean value of Mdot remain constant in %e\n",MDOT_ERR);        
         #ifdef _OPENMP                                                                  
         int time_sec = (int)(omp_get_wtime()-start);
         int hr  = time_sec/3600;
         int min = (time_sec%3600)/60;
         int sec = (time_sec%60)%60;
         printf("Expend %d hr : %d min : %d sec in the parallelized version using %d threads of %d available.\n",hr,min,sec,OMP_NUM,MAX_NUM_THREADS);
         #else                                                                     
         int time_sec = (int)((double)(clock()-start)/CLOCKS_PER_SEC);
         int hr  = time_sec/3600;
         int min = (time_sec%3600)/60;
         int sec = (time_sec%60)%60;
         printf("Expend %d hr : %d min : %d sec in serial version.\n",hr,min,sec);
         #endif                                                                    
         exit(EXIT_FAILURE);                                                       
      }    
      #endif
   }
}
