#include"main.h"
    
void Cons2Prim(double *u, double *q, gauge_ local_grid)
{
   int i, j, k;
   double D, E, S1=0.0, S2=0.0, S3=0.0;
   
   D  = q[DEN];
   E  = q[ENE];

#if DIM == 1
   S1 = q[MX1];
#elif DIM == 2
   S1 = q[MX1];
   S2 = q[MX2];
<<<<<<< HEAD
#elif DIM == 2
=======
#elif DIM == 3 || DIM == 4
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
   S1 = q[MX1];
   S2 = q[MX2];
   S3 = q[MX3];
#endif

   u[RHO] = D;
   #if EOS == IDEAL
   u[PRE] = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
   #elif EOS == DUST
   u[PRE] = 0.0;
   #endif
   u[VX1] = S1/D;
   u[VX2] = S2/D;
   u[VX3] = S3/D;
}
