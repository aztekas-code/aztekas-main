/**
 * @file eos.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Ultra-relativistic stiff Equation of state
 *
 */

#include"main.h"

void EoS(eos_ *eos, double *u, gauge_ local_grid)
{
   double rho, p;
   rho = u[EOS];
   p   = u[PRE];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   eos->h  = eos->e + p/rho;
   eos->cs = sqrt(K * p / (rho * eos->h));
#endif
}
