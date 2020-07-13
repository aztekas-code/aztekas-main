/**
 * @file eos.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Equation of state - Real relativistic gas (Ryu 2006) 
 * https://arxiv.org/pdf/astro-ph/0605550.pdf
 *
 */

#include"main.h"

void EoS(eos_ *eos, double *u, gauge_ *local_grid)
{
   double rho, p, theta_0;
   rho = u[RHO];
   p   = u[PRE];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   theta_0 = p / rho;
   eos->h = 2.0*((6.0*theta_0*theta_0 + 4.0*theta_0 + 1.0)/(3.0*theta_0 + 2.0));
   eos->cs = sqrt((theta_0*(3.0*theta_0 + 2.0)*(18.0*theta_0*theta_0 + 24.0*theta_0 + 5.0))/(3.0*(6.0*theta_0*theta_0 + 4.0*theta_0 + 1.0)*(9.0*theta_0*theta_0 + 12.0*theta_0 + 2.0)));
#endif
}
