/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:45:48
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

   Periodic(B);

   return 0;
}
