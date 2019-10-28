/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 18-10-2019 18:21:45
 * Created By :
 */

#include"main.h"

void Boundaries(double *B)
{
#if outflow_x1max == TRUE \
 || outflow_x1min == TRUE \
 || outflow_x2max == TRUE \
 || outflow_x2min == TRUE \
 || outflow_x3max == TRUE \
 || outflow_x3min == TRUE 

   Outflow(B);

#elif reflective_x1max == TRUE \
   || reflective_x1min == TRUE \
   || reflective_x2max == TRUE \
   || reflective_x2min == TRUE \
   || reflective_x3max == TRUE \
   || reflective_x3min == TRUE 

   Reflection(B);

#elif periodic_x1 == TRUE \
   || periodic_x2 == TRUE \
   || periodic_x3 == TRUE

   Periodic(B);

#endif
}
