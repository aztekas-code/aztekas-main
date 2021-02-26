/**
 * @file equation-solver.c
 *
 * @brief Define the type of system of equations to be solved.
 * Currently implemented:
 *    - HYDRO: Euler inviscid system of equations.
 *    - TOV: Tolman-Oppenheimer-Volkoff
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 26-02-2021 - 01:20:23
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 16-04-2020 - 17:58:21
 */
#include"main.h"

void Equation_System_Solver()
{
   int time_sec;
   int hr, min, sec;

#if TOV == TRUE

   ODE_Integration();

#endif

#if HYDRO == TRUE
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
      Print_Time_Values(&tprint,&dtprint,&itprint);

      //In here we set the integration method (MoL-RK and HRSC)
      Hyperbolic_Integration();
   
   #ifdef _OPENMP
      time_sec = (int)(omp_get_wtime() - start);
      hr       = time_sec/3600;
      min      = (time_sec%3600)/60;
      sec      = (time_sec%60)%60;
   #else
      time_sec = (int)((double)(clock()-start)/CLOCKS_PER_SEC);
      hr       = time_sec/3600;
      min      = (time_sec%3600)/60;
      sec      = (time_sec%60)%60;
   #endif
   #if PRINT_EVOLV == TRUE
      printf("Time = %e, dt = %e, Running time = %d hr : %d min : %d sec\r",\
            grid.time,dt,hr,min,sec);
      fflush(stdout);
   #endif
   }

   Print_Time_Values(&tprint,&dtprint,&itprint);
#endif
}
