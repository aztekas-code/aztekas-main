/*
 * aztekas boundaries module
 * Date of creation/modification: 26-09-19 11:33:21
 * author: Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
#if DIM == 1

   double r, r1,r2,rho1,rho2,m1, m2, mass_aux;

   // update mass vector 
   r1 = grid.X1[gc];
   rho1 = B(RHO,gc);
   mass[gc] = 4.*M_PI*r1*r1*r1*rho1/3. ;   
   for(int i = gc+1; i < Nx1-gc; i++)
   {
        r1 = grid.X1[i-1];
        r2 = grid.X1[i];        
        rho1 = B(RHO,i-1);
        rho2 = B(RHO,i);        
        m1 = 4.*M_PI*r1*r1*rho1 ;
        m2 = 4.*M_PI*r2*r2*rho2 ;
        mass[i] = mass[i-1] + 0.5*(m1+m2)*(r2-r1);
   }
   // accreted mass so far
   mass_acc = mass_tot - mass[Nx1-gc-1];
   if (mass_acc < 0.0) { mass_acc = 0.0 ;}
   //printf( "accreted mass : %e \n", mass_acc );
   // shift mass vector by mass_acc
    for(int i = gc-1; i < Nx1-gc; i++)
   {
        mass[i] = mass[i] + mass_acc;
   }  

   for(int i = 0; i < Nx1; i++)
   {
         B(PRE,i) = pressure_0;
   }

  
   for(int i = Nx1-gc; i <= Nx1; i++)
   {
         B(RHO,i) = density_atm;
         B(PRE,i) = pressure_atm;
         B(VX1,i) = 0.0;
   }

#elif DIM == 2

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0q;
            B(VX2,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0;
            B(PRE,i,j) = pressure_0;
            B(VX1,i,j) = velocity_0;
            B(VX2,i,j) = 0.0;
            B(VX3,i,j) = 0.0;
         }
      }
   }

#endif
}
