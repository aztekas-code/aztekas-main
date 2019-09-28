/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Rayleigh-Taylor
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:01:28
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int n, i, j, k, cell;

   Periodic(B);
   Reflection(B);

   return 0;
}
