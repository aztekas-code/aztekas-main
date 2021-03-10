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

   xMass[0] = 0.75;
   xMass[1] = 0.23;
   xMass[2] = 0.02;

   A[0] = 1.0;
   A[1] = 4.0;
   A[2] = 12.0;

   Z[0] = 1.0;
   Z[1] = 2.0;
   Z[2] = 6.0;

   rho = u[0];
   p   = u[1];
   e   = u[2];

   // Density and Pressure (prim2cons and fluxes and sources)
   if (e == 0.0)
      nad_eos_dp_(&rho,&p,xMass,A,Z,term_var);

   // Density and Energy (cons2prim)
   if (p == 0.0)
   {
      nad_eos_de_(&rho,&e,xMass,A,Z,term_var);
   }

   eos->p  = term_var[1];
   eos->e  = term_var[2];
   eos->cs = term_var[5];
}
