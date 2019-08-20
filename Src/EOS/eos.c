#include"main.h"

#if EOS == IDEAL

void EoS(eos_ *eos, double *u, grid_ local_grid)
{
   double rho, p;
   rho = u[0];
   p   = u[1];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   eos->h       = 1.0 + eos->e + p/rho;
   eos->cs = sqrt(K * p / (rho * eos->h));
#endif
}

#elif EOS == DUST

void EoS(eos_ *eos, double *u, grid_ local_grid)
{
   double rho, p;
   rho = u[0];
   p   = 0.0;

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   eos->h       = 1.0 + eos->e + p/rho;
   eos->cs = sqrt(K * p / (rho * eos->h));
#endif
}

#endif
