/**
 * @file eos.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Helmholtz equation of state
 *
 */

#include"main.h"

void EoS(eos_ *eos, double *u, gauge_ *local_grid)
{
   double xMass[3], A[3], Z[3];
   double term_var[6];
   double rho, p, e;

#ifndef HELMHOLTZ_COMP
   printf("\n");
   printf("HELMHOLTZ_COMP not defined\n");
   printf("Include in user_param.h one of the following options:\n");
   printf("   - #define HELMHOLTZ_COMP    DEFAULT -> H (75%), He (23%), C (2%)\n");
   printf("   - #define HELMHOLTZ_COMP    CO1     -> C (50%), O(50%)\n");
   printf("   - #define HELMHOLTZ_COMP    CO2     -> C (30%), O(70%)\n");
   printf("\n");
   exit(EXIT_FAILURE);
#endif

#if HELMHOLTZ_COMP == DEFAULT
   xMass[0] = 0.75;
   xMass[1] = 0.23;
   xMass[2] = 0.02;

   A[0] = 1.0;
   A[1] = 4.0;
   A[2] = 12.0;

   Z[0] = 1.0;
   Z[1] = 2.0;
   Z[2] = 6.0;
#elif HELMHOLTZ_COMP == CO1
   xMass[0] = 0.0;
   xMass[1] = 0.50;
   xMass[2] = 0.50;

   A[0] = 1.0;
   A[1] = 12.0;
   A[2] = 16.0;

   Z[0] = 1.0;
   Z[1] = 6.0;
   Z[2] = 8.0;
#elif HELMHOLTZ_COMP == CO2
   xMass[0] = 0.0;
   xMass[1] = 0.30;
   xMass[2] = 0.70;

   A[0] = 1.0;
   A[1] = 12.0;
   A[2] = 16.0;

   Z[0] = 1.0;
   Z[1] = 6.0;
   Z[2] = 8.0;
#endif

   rho = u[0];
   p   = u[1];
   e   = u[2];

   // Density and Pressure (prim2cons and fluxes and sources)
   if (e == 0.0)
      nad_eos_dp_(u,xMass,A,Z,term_var);

   // Density and Energy (cons2prim)
   if (p == 0.0)
   {
      nad_eos_de_(u,xMass,A,Z,term_var);
   }

   eos->rho  = term_var[0];
   eos->p    = term_var[1];
   eos->e    = term_var[2];
   eos->s    = term_var[3];
   eos->temp = term_var[4];
   eos->cs   = term_var[5];
}
