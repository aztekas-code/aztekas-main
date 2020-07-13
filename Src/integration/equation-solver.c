/**
 * @file equation-solver.c
 *
 * @brief Define the type of system of equations to be solved.
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 23-04-2020 - 16:16:40
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 16-04-2020 - 17:58:21
 */
#include"main.h"

void Equation_System_Solver()
{
#if HYPERBOLIC == TRUE
   while(grid.time <= tmax)
   {
      /**
       * Computes the time-step
       * File: /Src/integration/timestep.c
       */
      dt = TimeStep();

      /**
       * Real-Time computations using U
       */

      //We print the values
      Print_Values(&tprint,&dtprint,&itprint);

      //In here we set the integration method (MoL-RK and HRSC)
      Hyperbolic_Integration();

      printf("Time = %e, dt = %e\r",grid.time,dt);
      fflush(stdout);
   }

   Print_Values(&tprint,&dtprint,&itprint);
#endif

#if ODE == TRUE

   ODE_Integration();

#endif
}
