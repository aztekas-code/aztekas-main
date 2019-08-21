/* 
 *  aztekas boundaries module
 *  Date of creation: 17-05-2019 19:07:17
 *  author: Alejandro Aguayo-Ortiz 
 */
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"main.h"
#include"param.h"

int BOUNDARIES(double *B)
{
   int n, i, j, k, cell;
   double r;

   OUTFLOW(B);

   for(i = 0; i <= Nx1; i++)
   {   
      r_in = X1[gc]; 
      r = X1[i];
                                                                  
      if(r >= r_out)                                                    
      {                                                                
         B[c1(0,i)] = density_0;                                     
         B[c1(1,i)] = pressure_0;                                    
         B[c1(2,i)] = velocity_0;
      }                                                                
      else if(r <= r_in)                                               
      {                                                                
         B[c1(0,i)] = B[c1(0,gc)];
         B[c1(1,i)] = B[c1(1,gc)];
         B[c1(2,i)] = B[c1(2,gc)];
      }
   }

   return 0;
}
