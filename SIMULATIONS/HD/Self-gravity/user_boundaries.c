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

   polyK = (K-1.)/K*(1.- r_acc/Rad)/r_acc ;
   
   pressure_0 = polyK*pow(density_0,K);
   
   mass_aux = 0.0;
   // update mass vector 
   for(int i = gc-1; i < Nx1-gc; i++)
   {
        r1 = grid.X1[i];
        r2 = grid.X1[i+1];        
        rho1 = B(RHO,i);
        rho2 = B(RHO,i+1);        
        m1 = 4.*M_PI*r1*r1*rho1 ;
        m2 = 4.*M_PI*r2*r2*rho2 ;
        mass[i] = mass_aux + 0.5*(m1+m2);  
        mass_aux = mass[i];
   }
   // accreted mass so far
   mass_acc = mass_tot - mass[Nx1-gc-1];
   
   // shift mass vector by mass_acc
    for(int i = gc-1; i < Nx1-gc; i++)
   {
        mass[i] = mass[i] + mass_acc;
   }  

  
   for(int i = 0; i <= Nx1; i++)
   {

/*
      if(i <= gc)
      {

		  r = grid.X1[i];
	      B(RHO,i) = pow(1. - Rad/r*((r-1.)/(Rad-1.)),1./(K-1.));
    	  B(PRE,i) = polyK*pow(B(RHO,i),K);
    	  B(VX1,i) = 0.0;
	  }
	  

      if(i <= gc)
      {

	      B(RHO,i) = B(RHO,gc);
    	  B(PRE,i) = B(PRE,gc);
    	  B(VX1,i) = B(VX1,gc);
    	  if (B(VX1,gc) > 0) 
    	  {
    	    B(VX1,i) = 0.0;
    	  }
	  }
*/

      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_atm;
         B(PRE,i) = pressure_atm;
         B(VX1,i) = 0.0;
      }
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
