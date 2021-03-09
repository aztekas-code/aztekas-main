/*
 * aztekas initial module
 * Date of creation/modification: 27-09-19 11:26:00
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

mass = (double *)malloc((Nx1+1)*sizeof(double));

void Initial()
{
   int n, i, j, k, cell;
   
   double r;
   
   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;
   
   polyK = (K-1.)/K*(1.- r_acc/Rad)/r_acc/pow(density_0,K-1.) ;
   
#if DIM == 1 

   ///////////////////////////
   //-------Bondi-1D--------//
   ///////////////////////////
   
   double test, mass_aux;
   double density_atm, pressure_atm;
   
   int iRad;
   
   r=Rad;
   test = (Nx1-2*gc)*log(r-x1min+1)/log(x1max-x1min+1) + gc;
   
   iRad = (int)floor(test);
      
   //printf("%d %f \n",iRad,grid.X1[iRad]);
   
   for(i = 0; i <= iRad; i++)
   {
   	  r = grid.X1[i];
      U(RHO,i) = density_0*pow(1. - Rad/r*((r-r_acc)/(Rad-r_acc)),1./(K-1.));
      U(PRE,i) = polyK*pow(U(RHO,i),K);
      U(VX1,i) = 0.0;            
   }
   
   // define atmosphere density as 1000 times less 
   // than the minimum cloud's density 
   density_atm = 1.e-3*U(RHO,iRad);
   pressure_atm = polyK*pow(density_atm,K);

   for(i = iRad+1; i <= Nx1; i++)
   {
      r = grid.X1[i];
      U(RHO,i) = density_atm;
      U(PRE,i) = pressure_atm;
      U(VX1,i) = 0.0;
   }

   //initialize gas cloud mass vector
   mass_aux = 0.0;
   for(int i = gc-1; i < Nx1-gc; i++)
   {
        r1 = grid.X1[i];
        r2 = grid.X1[i+1];        
        rho1 = U(RHO,i);
        rho2 = U(RHO,i+1);        
        m1 = 4.*M_PI*r1*r1*rho1 ;
        m2 = 4.*M_PI*r2*r2*rho2 ;
        mass[i] = mass_aux + 0.5*(m1+m2);  
        mass_aux = mass[i];
   }
   //initial gas cloud's mass
   mass_tot = mass[Nx1-gc-1];
   mass_acc = 0.0;

#elif DIM == 2

   ///////////////////////////
   //-------Bondi-2D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = velocity_0;
         U(VX2,i,j) = 0.0;
      }
   }

#elif DIM == 4

   ///////////////////////////
   //------Bondi-2.5D-------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) = density_0;
         U(PRE,i,j) = pressure_0;
         U(VX1,i,j) = velocity_0;
         U(VX2,i,j) = 0.0;
         U(VX3,i,j) = 0.0;
      }
   }

#endif
}
