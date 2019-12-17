<<<<<<< HEAD
/*
 * File Name : user_boundaries.c
 * Description : aztekas boundaries module for Michel Stiff EoS problem
 * Creation Date : 05-12-2019
 * Last Modified : 05-12-2019 13:08:56
 * Created By : Alejandro Aguayo-Ortiz
 */

=======
>>>>>>> master
#include"main.h"

void User_Boundaries(double *B)
{
   int i, j, k, n, cell;
   double r, theta;
   double Delta, Sigma, rho2;
   double M, a;
   double rplus, rminus;

   Outflow(B);
#if POLAR == FALSE
   Reflection(B);
#elif POLAR == TRUE
   Periodic(B);
#endif

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(0,i) = density_0;
         B(1,i) = pressure_0;
         B(2,i) = velocity_0;
      }
   }

#elif DIM == 2

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            r     = grid.X1[i];
            theta = grid.X2[j];
            M     = Black_Hole_Mass;
            a     = Black_Hole_Spin;
         
            Delta  = r*r - 2.0*M*r + a*a;
            Sigma  = pow(r*r + a*a,2.0) - Delta*a*a*pow(sin(theta),2.0);
            rho2   = r*r + a*a*pow(cos(theta),2.0);
            rplus  = M + sqrt(M*M - a*a);
            rminus = M - sqrt(M*M - a*a);
         
            B(RHO,i,j) = sqrt(1 + ((2*M)/(rho2))*((r*(r + rplus) + 2*M*rplus)/(r-rminus)));
            B(PRE,i,j) = pow(B(RHO,i,j),K);
            B(VX1,i,j) = -(rplus*rplus + a*a)*sqrt(rho2/(Delta*Sigma));
            B(VX2,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         if(i >= Nx1-gc)
         {
            B(0,i,j) = density_0;
            B(1,i,j) = pressure_0;
            B(2,i,j) = velocity_0;
            B(3,i,j) = 0.0;
            B(4,i,j) = 0.0;
         }
      }
   }

#endif
}
