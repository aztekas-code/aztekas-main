/**
 * @file auxfunc.c
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Helpful functions for \a aztekas.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"
#include"param.h"

int MxV(double *M, double *V, double *L)
{
   int n, m;
   double res=0.0;

   for(m = 0; m < eq; m++)
   {
      for(n = 0; n < eq; n++)
      {
         res += M[m*(eq) + n]*V[n];
      }

      L[m] = res;
      res = 0.0;
   }

   return 0;
}

void RoundGen(double *num)
{
   double r;
   double bla;
   double decnum;
   int expnum;
   
   if(*num != 0.0e+00)
   {
      expnum = (int)floor(log(fabs(*num))/log(10.0));
      decnum = *num/pow(10,(double)expnum);
      decnum = roundf(decnum*1.0e+15)/1.0e+15;
      *num = decnum*pow(10,(double)expnum);
   }
}
