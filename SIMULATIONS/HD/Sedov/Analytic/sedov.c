/**
 * @mainpage SEDOV: analytic solution
 *
 * This program solves the Sedov-Taylor problem in a one dimensional case
 * for a ideal gas with a polytropic equation of state \f$ p \propto \rho^\gamma \f$.
 *
 * The solution follows the treatment presented in Fluid Mechanics of
 * Landau \& Lifschitz, Second Edition, Volume 6 of Course of Theoretical
 * Physics (1987).
 */

/**
* @file sedov.c
* @author Alejandro Aguayo-Ortiz
* @date 30 Mar 2019
* @brief Main program
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>

/**
* @brief Main function
*/
int main(int argc, char* argv[])
{

   /** @brief Define dummie variables */
   int i;
   /** @brief Define input variables */
   int N = 512;
   int coord = 3;
   double *Xi;
   double g = 5.0/3.0;

   Xi = (double *)malloc((N+1) * sizeof(double));

   /* @brief Define the physical radius \f$ r \f$*/
   for(i = 0; i <= N; i++)
   {
      Xi[i] = (double)i / (double)N;
   }

   // Define useful variables
   double nu1, nu2, nu3, nu4, nu5;

   nu1 = -(13 * pow(g,2.0) - 7 * g + 12)/((3 * g - 1) * (2 * g + 1));
   nu2 = (5 * (g - 1))/(2 * g + 1);
   nu3 = 3/(2 * g + 1);
   nu4 = - nu1 / (2 - g);
   nu5 = - 2 /(2 - g);

   /** @brief Define adimentional variables V, G and Z */
   double *V, *G, *Z;
   V = (double *)malloc((N + 1) * sizeof(double));
   G = (double *)malloc((N + 1) * sizeof(double));
   Z = (double *)malloc((N + 1) * sizeof(double));

   /** Root finding */
   double xi, v, v0, f, df;
   double term1, term2, term3;
   double dterm1, dterm2, dterm3;

   /** Useful variables */
   double gp1 = g + 1;
   double gm1 = g - 1;
   double gpm = gp1 / gm1;
   double gp7 = gp1 / (7 - g);
   double g13 = 3 * g - 1;

   for(i = 0; i <= N; i++)
   {
      xi = Xi[i];
      v0 = 2.0 / gp1;
      f  = 1.0;
      df = 1.0;
      
      while(fabs(f) > 0.0001)
      {
         term1 = 4.0 / (pow(gp1 * v0,2.0));
         term2 = pow(gp7,nu1) * pow(5 - g13 * v0,nu1);
         term3 = pow(gpm,nu2) * pow(g * v0 - 1,nu2);

         dterm1 = - 8.0 / (pow(gp1,2.0) * pow(v0,3.0));
         dterm2 = pow(gp7,nu1) * nu1 * (-g13) * pow(5 - g13 * v0,nu1 - 1.0);
         dterm3 = pow(gpm,nu2) * nu2 * g * pow(g * v0 - 1,nu2 - 1.0);

         f  = term1 * term2 * term3 - pow(xi,5.0);
         df = dterm1 *  term2 *  term3 + \
               term1 * dterm2 *  term3 + \
               term1 *  term2 * dterm3;

         v = v0 - f/df;

         if(v <= 1.0 / g)
         {
            v = v0 - 0.5*f/df;
         }

         v0 = v;

      }

      V[i] = v0;
      Z[i] = g * gm1 * (1 - V[i]) * pow(V[i],2.0) / (2 * (g * V[i] - 1));
      G[i] = gpm * pow(gpm * (g * V[i] - 1),nu3) * pow(gp7 * (5 - g13 * V[i]),nu4) * pow(gpm * (1 - V[i]),nu5);
   }

   /** Computing beta */
   double integral = 0;
   double beta;

   for(i = 0; i <= N; i++)
   {
      integral = integral + G[i] * (0.5 * pow(V[i],2.0) + Z[i] / (g * gm1)) * pow(Xi[i],4.0) / N;
   }

   beta = pow(25 / (16 * M_PI * integral),1.0/5.0);

   double time, Rmax, Rmin, E, rho1, pre1, vel1, R;
   /** Recover physical quantities*/
   time  = 0.1;
   Rmax  = 8.0;
   Rmin  = 0.0;
   E     = 1.25e+05;
   rho1  = 1.0;
   pre1  = 1.0e-05;
   vel1  = 0.0;

   /** Shock position*/
   R = beta*pow(E * pow(time,2.0)/rho1,0.2);

   double *r;
   r = (double *)malloc((N+1) * sizeof(double));

   /* @brief Define the physical radius \f$ r \f$*/
   for(i = 0; i <= N; i++)
   {
      r[i]  = (double)i*(Rmax-Rmin)/((double)N);
      Xi[i] = r[i]/R;
   }

   for(i = 0; i <= N; i++)
   {
      xi = Xi[i];
      v0 = 2.0 / gp1;
      f  = 1.0;
      df = 1.0;
      
      while(fabs(f) > 0.0001)
      {
         term1 = 4.0 / (pow(gp1 * v0,2.0));
         term2 = pow(gp7,nu1) * pow(5 - g13 * v0,nu1);
         term3 = pow(gpm,nu2) * pow(g * v0 - 1,nu2);

         dterm1 = - 8.0 / (pow(gp1,2.0) * pow(v0,3.0));
         dterm2 = pow(gp7,nu1) * nu1 * (-g13) * pow(5 - g13 * v0,nu1 - 1.0);
         dterm3 = pow(gpm,nu2) * nu2 * g * pow(g * v0 - 1,nu2 - 1.0);

         f  = term1 * term2 * term3 - pow(xi,5.0);
         df = dterm1 *  term2 *  term3 + \
               term1 * dterm2 *  term3 + \
               term1 *  term2 * dterm3;

         v = v0 - f/df;

         if(v <= 1.0 / g)
         {
            v = v0 - 0.5*f/df;
         }

         v0 = v;

      }

      V[i] = v0;
      Z[i] = g * gm1 * (1 - V[i]) * pow(V[i],2.0) / (2 * (g * V[i] - 1));
      G[i] = gpm * pow(gpm * (g * V[i] - 1),nu3) * pow(gp7 * (5 - g13 * V[i]),nu4) * pow(gpm * (1 - V[i]),nu5);
   }

   double *rho, *vr, *c, *p;
   rho = (double *)malloc((N + 1) * sizeof(double));
   vr  = (double *)malloc((N + 1) * sizeof(double));
   c   = (double *)malloc((N + 1) * sizeof(double));
   p   = (double *)malloc((N + 1) * sizeof(double));

   for(i = 0; i <= N; i++)
   {
      rho[i] = rho1*G[i];
      c[i]   = sqrt(4.0*r[i]*r[i]*Z[i]/(25.0*time*time));
      p[i]   = rho[i]*c[i]*c[i]/g;
      vr[i]  = 2.0 * r[i] * V[i] / (5.0 * time);    
      if(Xi[i] > 1.0)
      {
         rho[i] = rho1;
         p[i]   = pre1;
         vr[i]  = vel1;
      }
      printf("%e %e %e %e \n",r[i], rho[i], p[i], vr[i]);
   }

   return 0;
}
