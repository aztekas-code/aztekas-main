/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 28-10-2019 17:31:22
 * Created By :
 */

#include"main.h"

void User_Boundaries(double *B)
{
#if DIM == 2
   
   // Linear extrapolation. Important for the diagonal shock tube
   for(int n = 0; n < eq; n++)
   {
      for(int i = gc; i <= Nx1-gc; i++)
      {
         for(int j = gc; j <= Nx2-gc; j++)
         {
            if(i > Nx1-gc)
            {
               B(n,i,j) = B(n,i-2,j) + ((grid.X1[i] - grid.X1[i-2])/(grid.X1[i-1] - grid.X1[i-2]))*(B(n,i-1,j) - B(n,i-2,j)); 
            }
            if(j > Nx2-gc)
            {
               B(n,i,j) = B(n,i,j-2) + ((grid.X2[j] - grid.X2[j-2])/(grid.X2[j-1] - grid.X2[j-2]))*(B(n,i,j-1) - B(n,i,j-2)); 
            }
         }
      }

      for(int i = Nx1-gc; i >= gc; i--)
      {
         for(int j = Nx2-gc; j >= gc; j--)
         {
            if(i < gc)
            {
               B(n,i,j) = B(n,i+2,j) + ((grid.X1[i] - grid.X1[i+2])/(grid.X1[i+1] - grid.X1[i+2]))*(B(n,i+1,j) - B(n,i+2,j)); 
            }
            if(j < gc)
            {
               B(n,i,j) = B(n,i,j+2) + ((grid.X2[j] - grid.X2[j+2])/(grid.X2[j+1] - grid.X2[j+2]))*(B(n,i,j+1) - B(n,i,j+2)); 
            }
         }
      }
   }
#endif
}
