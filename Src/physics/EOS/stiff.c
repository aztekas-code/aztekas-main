/**
 * @file eos.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Equation of state
 *
 */

#include"main.h"

void EoS(eos_ *eos, double *u, gauge_ local_grid)
{
   double rho, p;
   rho = u[0];
   p   = u[1];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   eos->h  = eos->e + p/rho;
   eos->cs = sqrt(K * p / (rho * eos->h));
#endif
}
