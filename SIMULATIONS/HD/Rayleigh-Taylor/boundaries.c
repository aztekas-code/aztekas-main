/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 12:50:45
 *  author: Alejandro Aguayo Ortiz 
 */
#include"main.h"

int Boundaries(double *B)
{
   int n, i, j, k, cell;

   Periodic(B);
   Reflection(B);

   return 0;
}
