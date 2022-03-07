/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 17:40:31
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
      {
         if(i <= gc)
         {
            U(RHO,i,j) = 0.01;
            U(PRE,i,j) = 0.0015;
            U(VX1,i,j) = 1.0;
            U(VX2,i,j) = 0.0;
         }
      }
   }
}
