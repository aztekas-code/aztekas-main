/* 
 *  aztekas initial module
 *  Date of creation: 02-01-2019 18:32:00
 *  author: Alejandro Aguayo Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include"main.h"
#include"param.h"

void INITIAL(double *dtprint)
{
   int n, i, j, k, cell;

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(X1[i] <= r_jet && X2[j] <= z_jet)
         {
            U[c2(0,i,j)] = n_jet;
            U[c2(1,i,j)] = p_jet;
            U[c2(2,i,j)] = vx1_jet;
            U[c2(3,i,j)] = vx2_jet;
         }

         U[c2(0,i,j)] = n_atm;
         U[c2(1,i,j)] = p_atm;
         U[c2(2,i,j)] = vx1_atm;
         U[c2(3,i,j)] = vx2_atm;
      }
   }
}
