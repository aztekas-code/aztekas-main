/*
 * aztekas initial module
 * Date of creation: 26-09-2019 00:09:12
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

int Boundaries(double *B)
{
   int i, j, k, n, cell;

   Outflow(B);
   Reflection(B);

   return 0;
}
