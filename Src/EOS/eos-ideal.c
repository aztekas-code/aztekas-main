#include"main.h"

void EoS_Ideal(eos_ *eos, double *u, double *x)
{
   double rho, p, h;
   rho = u[0];
   p   = u[1];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   h       = 1.0 + eos->e + p/rho;
   eos->cs = sqrt(K * p / (rho * h));
#endif
}
