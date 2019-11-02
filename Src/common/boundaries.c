/**
 * @file boundaries.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Standard boundary conditions. Outflow, Periodic and Reflection.
 */

#include"main.h"

/**
 * Call the default boundary conditions.
 *
 * Input/Output: double *B
 */
void Boundaries(double *B)
{
#if reflective_x1max == TRUE \
 || reflective_x1min == TRUE \
 || reflective_x2max == TRUE \
 || reflective_x2min == TRUE \
 || reflective_x3max == TRUE \
 || reflective_x3min == TRUE 

   Reflection(B);

#endif

#if periodic_x1 == TRUE \
 || periodic_x2 == TRUE \
 || periodic_x3 == TRUE

   Periodic(B);

#endif

#if outflow_x1max == TRUE \
 || outflow_x1min == TRUE \
 || outflow_x2max == TRUE \
 || outflow_x2min == TRUE \
 || outflow_x3max == TRUE \
 || outflow_x3min == TRUE 

   Outflow(B);

#endif

   User_Boundaries(B);
}
