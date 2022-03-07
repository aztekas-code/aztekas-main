/*
 * File Name : initial.c
 * Description : aztekas initial module for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 07-10-2019 21:23:58
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   //////////////////////////////////
   // Kelvin-Helmholtz Instability //
   //////////////////////////////////
   for(int i = 0; i <= Nx1; i++)                                                    
   {                                                                            
      for(int j = 0; j <= Nx2; j++)                                                 
      {                                                                         
         double r = sqrt(pow(grid.X1[i] + 0.5,2.0) + pow(grid.X2[j],2.0));                
                                                                                
         U(PRE,i,j) = 0.0015;                                                 
         U(VX2,i,j) = 0.0;                                                    
                                                                                
         if(r <= r_in)                                                          
         {                                                                      
            U(RHO,i,j) = 100.0;                                               
            U(VX1,i,j) = 0.0;                                                 
         }                                                                      
         else
         {                                                                      
            U(RHO,i,j) = 0.01;                                                
            U(VX1,i,j) = 1.0;                                                 
         }                                                                      
      }                                                                         
   }  
   //////////////////////////////////
}
