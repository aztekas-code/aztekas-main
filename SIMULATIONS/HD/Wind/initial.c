/*
 * aztekas initial module
 * Date of creation/modification: 27-09-19 11:26:00
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   double u[3];
   eos_ eos;
   gauge_ local_grid;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   ///////////////////////////
   //---------Wind----------//
   ///////////////////////////

#ifdef HELMHOLTZ_COMP
   u[0] = density_0;
   u[1] = temperature_0;
   u[2] = 0.0;

   EoS_DT(&eos,u);

   pressure_0 = eos.p;
#endif

   G = 6.67e-08; 
   M = 2.0e+33;

   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) =  density_0;
         U(PRE,i,j) =  pressure_0;//pow(U(RHO,i,j),K)/K;
         U(VX1,i,j) =  velocity_0*cos(grid.X2[j]);
         U(VX2,i,j) = -velocity_0*sin(grid.X2[j]);
      }
   }
}
