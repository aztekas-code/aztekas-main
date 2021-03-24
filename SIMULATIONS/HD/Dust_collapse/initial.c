/*
 * aztekas initial module
 * Date of creation/modification: 27-09-19 11:26:00
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;
   
   double r;
   
   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;
   

#if DIM == 1 

   ///////////////////////////
   //-------Bondi-1D--------//
   ///////////////////////////
   
   double test, mass_aux;
   double r1,r2,rho1,rho2,m1,m2;
   
   int iRad;
   
   r=Rad;
   test = (Nx1-2*gc)*log(r-x1min+1)/log(x1max-x1min+1) + gc;
   
   iRad = (int)floor(test);
      
   //printf("%d %f \n",iRad,grid.X1[iRad]);
   
   for(i = 0; i <= iRad; i++)
   {
   	  r = grid.X1[i];
   	  U(RHO,i) = density_0;
   	  U(PRE,i) = pressure_0;
      U(VX1,i) = 0.0;
   }
   
   // define atmosphere density as 1000 times less 
   // than the minimum cloud's density 
   //printf( "%e %e \n", grid.X1[iRad], U(RHO,iRad) );   
   density_atm = 1e-3*density_0;
   pressure_atm = pressure_0;

   for(i = iRad+1; i <= Nx1; i++)
   {
      r = grid.X1[i];
      U(RHO,i) = density_atm;
      U(PRE,i) = pressure_atm;
      U(VX1,i) = 0.0;
   }

   //initialize gas cloud mass vector
   r1 = grid.X1[gc];
   mass[gc] = 4.*M_PI*r1*r1*r1*density_0/3. ;
   //printf( "%e, %e \n", r1, mass[gc] );
   
   for(int i = gc+1; i < Nx1-gc; i++)
   {
        r1 = grid.X1[i-1];
        r2 = grid.X1[i];        
        rho1 = U(RHO,i-1);
        rho2 = U(RHO,i);        
        m1 = 4.*M_PI*r1*r1*rho1 ;
        m2 = 4.*M_PI*r2*r2*rho2 ;
        mass[i] = mass[i-1] + 0.5*(m1+m2)*(r2-r1);  
        //printf( "%e %e %e \n", r2, mass[i], m2*r2/3. );
        //printf( "%e \n", r1 );
   }
   //initial gas cloud's mass
   mass_tot = mass[Nx1-gc-1];
   //printf( "%e \n", mass_tot );
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
