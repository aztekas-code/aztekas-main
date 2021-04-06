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
   double term_var[8];
   double rho, p, e;

#ifndef HELMHOLTZ_COMP
   printf("\n");
   printf("HELMHOLTZ_COMP not defined\n");
   printf("Include in user_param.h one of the following options:\n");
   printf("   - #define HELMHOLTZ_COMP    DEFAULT -> H (75%), He (23%), C (2%)\n");
   printf("   - #define HELMHOLTZ_COMP    CO1     -> C (50%), O(50%)\n");
   printf("   - #define HELMHOLTZ_COMP    CO2     -> C (30%), O(70%)\n");
   printf("\n");
   printf("Use function EoS_DT(eos_ *eos, double *u) in order to obtain\n");
   printf("the thermodynamical variables from Density = u[0] and Temperature = u[1]\n");
   printf("   - eos.rho\n");
   printf("   - eos.p\n");
   printf("   - eos.e\n");
   printf("   - eos.s\n");
   printf("   - eos.temp\n");
   printf("   - eos.cs\n");
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
   
   //u[0] = roundf(u[0] * 10)/ 10;
   //u[1] = roundf(u[1] * 10)/ 10;
   //u[2] = roundf(u[2] * 10)/ 10;
   
   //printf("%e\n",u[1]);

   u[0] = u[0]*dens_units;
   u[1] = u[1]*dens_units*vel_units*vel_units;
   u[2] = u[2]*vel_units*vel_units;

   // Density and Pressure (prim2cons and fluxes and sources)
   if (u[2] == 0.0)
      nad_eos_dp_(u,xMass,A,Z,term_var);

   // Density and Energy (cons2prim)
   if (u[1] == 0.0)
   {
      nad_eos_de_(u,xMass,A,Z,term_var);
   }

   eos->rho  = term_var[0]/(dens_units);
   eos->p    = term_var[1]/(dens_units*vel_units*vel_units);
   eos->e    = term_var[2]/(vel_units*vel_units);
   eos->s    = term_var[3]/(vel_units*vel_units/temp_units);
   eos->temp = term_var[4]/(temp_units);
   eos->cs   = term_var[5]/(vel_units);
   eos->dpt  = term_var[6]/(dens_units*vel_units*vel_units/temp_units);
   eos->det  = term_var[7]/(vel_units*vel_units/temp_units);
}

void EoS_DT(eos_ *eos, double *u)
{
   double xMass[3], A[3], Z[3];
   double term_var[8];
   double rho, p, e, temp;

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

   u[0] = u[0]*dens_units;
   u[1] = u[1]*temp_units;

   nad_eos_dt_(u,xMass,A,Z,term_var);

   eos->rho  = term_var[0]/(dens_units);
   eos->p    = term_var[1]/(dens_units*vel_units*vel_units);
   eos->e    = term_var[2]/(vel_units*vel_units);
   eos->s    = term_var[3]/(vel_units*vel_units/temp_units);
   eos->temp = term_var[4]/(temp_units);
   eos->cs   = term_var[5]/(vel_units);
   eos->dpt  = term_var[6]/(dens_units*vel_units*vel_units/temp_units);
   eos->det  = term_var[7]/(vel_units*vel_units/temp_units);
}
