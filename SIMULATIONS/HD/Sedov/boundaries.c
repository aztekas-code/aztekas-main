/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Sedov Blast Wave
 * Creation Date : 26-09-2019
 * Last Modified : 26-09-2019 23:56:42
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

   Outflow(B);
   Reflection(B);

   return 0;
}
