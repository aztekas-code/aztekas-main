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

//   for(n = 0; n < eq; n++)
//   {
//      for(i = 0; i <= Nx1; i++)
//      {
//         for(j = 0; j <= Nx2; j++)
//         {
//            RoundGen(&B(n,i,j));
//         }
//      }
//   }

   return 0;
}
