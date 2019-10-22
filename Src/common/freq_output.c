/*
 * File Name : frec_output.c
 * Description :
 * Creation Date : 22-10-2019
 * Last Modified : 22-10-2019 15:26:33
 * Created By :
 */

#include"main.h"

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
