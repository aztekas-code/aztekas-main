/**
 * @file freq_output.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Define the type of frequency on the output: number of 
 * files (numfile) or time interval (timefile)
 */

#include"main.h"

/**
 * This function gives the kind of frequency output choose by the user. 
 * Wheter it is by number of files or time interval.
 *
 * Output: double *dtprint
 */

void Frequency_Output(double *dtprint)
{
   if(timefile >= 0.0)
   {
      *dtprint = timefile;
   }

   if(numfile > 0)
   {
      *dtprint = (tmax - grid.time)/numfile;
   }
}
