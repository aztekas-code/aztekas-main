/**
 * @file integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

//Do not erase any of these libraries//
#include"main.h"

void Integration()
{
   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
#if DIM == 1 
   
   Prim2Cons_All(Q,U);

   RK1D(U,Q,Q1,Q2,1);
   Cons2Prim(U,Q1);
   Boundaries(U);

   RK1D(U,Q,Q1,Q2,2);
   Cons2Prim(U,Q2);
   Boundaries(U);
   
#elif DIM == 2 || DIM == 4

   Prim2Cons_All(Q,U);

   RK2D(U,Q,Q1,Q2,1);
   Cons2Prim(U,Q1);
   Boundaries(U);

   U0 = U;

   RK2D(U,Q,Q1,Q2,2);
   Cons2Prim(U,Q2);
   Boundaries(U);

   U0 = U;
   
#elif DIM == 3 
   
   Prim2Cons_All(Q,U);
   RK3D(U,Q,Q1,Q2,1);
   Cons2Prim(U1,Q1);
   Boundaries(U1);
   RK3D(U1,Q,Q1,Q2,2);
   Cons2Prim(U,Q2);
   Boundaries(U);
   
#endif
}
