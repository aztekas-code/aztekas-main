/* 
 *  aztekas boundaries module
 *  Date of creation: 02-01-2019 18:32:00
 *  author: Alejandro Aguayo Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

int BOUNDARIES(double *B)
{
   int n, i, j, k, cell;

   OUTFLOW(B);
   REFLECTIVE(B);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(X1[i] <= r_jet && X2[j] <= z_jet)
         {
            B[c2(0,i,j)] = n_jet;
            B[c2(1,i,j)] = p_jet;
            B[c2(2,i,j)] = vx1_jet;
            B[c2(3,i,j)] = vx2_jet;
         }

         B[c2(0,i,j)] = n_atm;
         B[c2(1,i,j)] = p_atm;
         B[c2(2,i,j)] = vx1_atm;
         B[c2(3,i,j)] = vx2_atm;
      }
   }

   return 0;
}
