/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 00:05:21
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

   Periodic(B);

   return 0;
}
