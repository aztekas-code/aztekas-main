/**
 * @file /fluxes/limiters.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Define the type of reconstruction.
 * GODUNOV, MINMOD, MC, SUPERBEE and WENO5 implemented.
 */

#include"main.h"

double Limiter(double A, double B, int reconst)
{
   double sig;

   sig = 0.0;

   switch(reconst)
   {
      case GODUNOV:
         sig = 0.0;
      break;

      case MINMOD:
         sig = Minmod(A,B);
      break;

      case MC:
         sig = Mc(A,B);
      break;

      case SUPERBEE:
         sig = Superbee(A,B);
      break;
   }

   return sig;
}

double Maxmod(double A, double B)
{
   double sig;

   if(A*B < 0)
   {
      sig = 0;
   }
   else if(A*B >= 0)
   {
      if(fabs(A) > fabs(B))
      {
         sig = A;
      }
      else if(fabs(A) < fabs(B))
      {
         sig = B;
      }
   }

   return sig;
}

double Minmod(double A, double B)
{
   double sig;

   if(A*B <= 0)
   {
      sig = 0;
   }
   else if(A*B > 0)
   {
      if(fabs(A) < fabs(B))
      {
         sig = A;
      }
      else if(fabs(A) > fabs(B))
      {
         sig = B;
      }
   }

   return sig;
}

double Mc(double A, double B)
{
   double sig;
   double C;

   C = (A + B) / 2.0;

   if(A*B <= 0)
   {
      sig = 0.0;
   }
   else if(A*B > 0)
   {
      if(fabs(A) < fabs(B) && 2.0*fabs(A) < fabs(C))
      {
         sig = 2.0*A;
      }
      else if(fabs(A) >= fabs(B) && 2.0*fabs(B) < fabs(C))
      {
         sig = 2.0*B;
      }
      else if(fabs(C) <= 2.0*fabs(A) && fabs(C) <= 2.0*fabs(B))
      {
         sig= C;
      }
   }

   return sig;
}

double Superbee(double A, double B)
{
   double sig1;
   double sig2;
   double sig;

   sig1 = Minmod(2*A,B);
   sig2 = Minmod(A,2*B);

   sig = Maxmod(sig1,sig2);
   return sig;
}

double Weno5(double v1, double v2, double v3, double v4, double v5)
{
   double c13d12 = 1.3e+01 / 1.2e+01;
   double c1d4   = 1.0e+00 / 4.0e+00;
   double c11d6  = 1.1e+01 / 6.0e+00;
   double c2d6  = 2.0e+00 / 6.0e+00;
   double c7d6  = 7.0e+00 / 6.0e+00;
   double c5d6  = 5.0e+00 / 6.0e+00;
   double c1d6  = 1.0e+00 / 6.0e+00;
   double c1d3  = 1.0e+00 / 3.0e+00;

   double fs1, fs2, fs3;
   double bs1, bs2, bs3;

   double ep = 1.0e-40;

   double P1, P2, P3;
   double s1, s2, s3;
   double a1, a2, a3;

   double asum;
   double g1, g2, g3;
   double d1, d2, d3;

   double om1, om2, om3;
   double R;

   d1 = 1.0e-01;
   d2 = 6.0e-01;
   d3 = 3.0e-01;

   P1 =  c1d3*v1 - c7d6*v2 + c11d6*v3;
   P2 = -c1d6*v2 + c5d6*v3 +  c1d3*v4;
   P3 =  c1d3*v3 + c5d6*v4 -  c1d6*v5;

   fs1 = (v1 - 2.0e+00 * v2 + v3);
   fs2 = (v2 - 2.0e+00 * v3 + v4);
   fs3 = (v3 - 2.0e+00 * v4 + v5);

   bs1 = (v1 - 4.0e+00*v2 + 3.0e+00*v3);
   bs2 = (v2 - v4);
   bs3 = (3.0e+00*v3 - 4.0e+00*v4 + v5);

   s1 = c13d12*fs1*fs1 + c1d4*bs1*bs1;
   s2 = c13d12*fs2*fs2 + c1d4*bs2*bs2;
   s3 = c13d12*fs3*fs3 + c1d4*bs3*bs3;

   a1 = d1 * (1.0e+00 / pow(s1 + ep,2.0e+00));
   a2 = d2 * (1.0e+00 / pow(s2 + ep,2.0e+00));
   a3 = d3 * (1.0e+00 / pow(s3 + ep,2.0e+00));

   asum = a1 + a2 + a3;

   om1 = a1 / asum;
   om2 = a2 / asum;
   om3 = a3 / asum;

   R = om1*P1 + om2*P2 + om3*P3;

   return R;
}
