/**
 * @file /integration/HYPERBOLIC/integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

#include"main.h"

void Hyperbolic_Integration()
{
   rk_order = 2;

   /**
    * Convert all primitive vector U to conservative Q0
    */
   Prim2Cons_All(Q0,U);
   
   for(int order = 1; order <= rk_order; order++)
   {
      Primitive_Reconstruction();
      Method_of_Lines(order);
      Cons2Prim(U,Q);
      Boundaries(U);
      U0 = U;
   }

   /**
    * Increase time by dt
    */
   grid.time = grid.time + dt;
}
