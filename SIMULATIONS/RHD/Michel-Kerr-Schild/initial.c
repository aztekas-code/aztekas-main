/*
 * aztekas initial module
 * Date of creation/modification: 26-09-19 11:56:23
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   int n, i, j, k, cell;
   double r, theta;
   double Delta, Sigma, rho2;
   double M, a;
   double rplus, rminus;
   double alpha, betar, grr, grp, gpp;
   double Ut, Ur, Up;
   double Vr, Vp, vr, vp;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   ////////////////////////////
   //-------Michel-1D--------//
   ////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      U(RHO,i) = density_0;
      U(PRE,i) = pressure_0;
      U(VX1,i) = velocity_0;
   }

#elif DIM == 2

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

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

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         r     = grid.X1[Nx1-gc];
         theta = grid.X2[Nx2-gc];
         M     = Black_Hole_Mass;
         a     = Black_Hole_Spin;
         
         Delta  = r*r - 2.0*M*r + a*a;
         Sigma  = pow(r*r + a*a,2.0) - Delta*a*a*pow(sin(theta),2.0);
         rho2   = r*r + a*a*pow(cos(theta),2.0);
         rplus  = M + sqrt(M*M - a*a);
         rminus = M - sqrt(M*M - a*a);

         alpha = pow(1.0 + 2.0*M*r/rho2,-0.5);
         betar = (2.0*M*r/rho2)*pow(1.0 + 2.0*M*r/rho2,-1.0);
         grr   = 1.0 + 2.0*M*r/rho2;
         grp   = -a*(1.0 + 2.0*M*r/rho2)*pow(sin(theta),2.0);
         gpp   = pow(sin(theta),2.0)*(rho2 + a*a*(1.0 + 2.0*M*r/rho2)*pow(sin(theta),2.0));

         Ut = 1.0 + 2.0*M*r*(r+rplus)/(rho2*(r-rminus));
         Ur = -2.0*M*rplus/rho2;
         Up = 2.0*M*a/(rho2*(r-rminus));

         Vr = Ur/Ut;
         Vp = Up/Ut;

         vr = Vr/alpha + betar/alpha;
         vp = Vp/alpha;

         U(RHO,i,j) = sqrt(1.0 + ((2.0*M)/(rho2))*((r*(r + rplus) + 2.0*M*rplus)/(r-rminus)));
         U(PRE,i,j) = pow(U(RHO,i,j),K);
         U(VX1,i,j) = grr*vr + grp*vp;
         U(VX2,i,j) = 0.0;
         U(VX3,i,j) = grp*vr + gpp*vp;
      }
   }

#endif
}
