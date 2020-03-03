#include"main.h"

void Initial()
{
   double r;
   double t;
   double vr;
   double alpha, beta;
   double rho2, Delta, Sigma;
   double rin = r_in;
   double rin2 = rin*rin;
   double q, qin, C;
   double term0, term1, term2;
   double term3, term5, term6;
   double term4, term4_1, term4_2;
   double M = Black_Hole_Mass;
   double a = Black_Hole_Spin;
   double lnh, h, e;
   double v_phi;
   double l;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   ////////////////////////////
   //-------Michel-1D--------//
   ////////////////////////////

   for(int i = 0; i <= Nx1; i++)
   {
      r = grid.X1[i];

      U(0,i) = rho(r);
      U(1,i) = pre(r);
      U(2,i) = velocity_0;
   }

#elif DIM == 2

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         r = grid.X1[i];

         U(0,i,j) = rho(r);
         U(1,i,j) = pre(r);
         U(2,i,j) = velocity_0;
         U(3,i,j) = 0.0;
      }
   }

#elif DIM == 4

   //////////////////////////////
   //-------Michel-2.5D--------//
   //////////////////////////////

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         r    = grid.X1[i];
         t    = grid.X2[j];
         l    = l_0;

         alpha = sqrt(r/(r + 2.0*M));
         beta  = (2.0*M/r);

         // In FM76, Delta = Delta, Sigma = rho2, A = Sigma
         q   = sqrt(1.0 + (4.0*l*l*(1.0 - 2.0*M/r))/pow(r*sin(t),2.0));
         qin = sqrt(1.0 + (4.0*l*l*(1.0 - 2.0*M/rin))/pow(rin,2.0));
         C   = 0.5*log((1.0 + qin)/(1.0 - 2.0*M/rin)) - 0.5*qin;

         lnh   = 0.5*log((1.0 + q)/(1.0 - 2.0*M/r)) - 0.5*q - C;
         h     = exp(lnh);
         e     = (h - 1.0)/K;
         K_pol = (K - 1.0)*pow(0.01,2.0)/((K - 1.0 - pow(0.01,2.0))*K);

         v_phi = (q - 1.0)*pow(r*sin(t),2.0)/(2.0*alpha*l);

         if(v_phi < 1.0e-05)
         {
            v_phi = 0.0;
         }

         U(VX3,i,j) = 0.0;

         if(r >= rin && e > 0.0)
         {
            U(RHO,i,j) = pow((K - 1.0)*e/K_pol,1.0/(K - 1.0));
            U(VX3,i,j) = v_phi;
         }
         else if(r > rin && e <= 0.0)
         {
            U(RHO,i,j) = 1.0e-04*pow(r/rin,-1.5);
         }
         else if(r < rin)
         {
            U(RHO,i,j) = 1.0e-04*pow(r/rin,-1.5);
         }

         U(PRE,i,j) = K_pol*pow(U(RHO,i,j),K);
         U(VX1,i,j) = beta/alpha;
         U(VX2,i,j) = 0.0;
      }
   }

#endif
}
