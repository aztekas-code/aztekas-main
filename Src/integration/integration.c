/**
 * @file /integration/integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

#include"main.h"

void Integration()
{
   int rk_order = 2;

   /**
    * Increase time by dt
    */
   grid.time = grid.time + dt;

   /**
    * Convert all primitive vector U to conservative Q0
    */
   Prim2Cons_All(Q0,U);
   
   for(int order = 1; order <= rk_order; order++)
   {
      Primitive_Reconstruction();
      Runge_Kutta(order);
      Cons2Prim(U,Q);
      Boundaries(U);
      U0 = U;
   }
}
