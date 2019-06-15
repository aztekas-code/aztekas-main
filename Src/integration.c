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

int Integration()
{
   int n, i, j, k;
   double uu[eq+1];
   double qq[eq+1];

   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
#if DIM == 1 
   
   Prim2Cons_All(Q,U);

   RK1D(U,Q,Q1,Q2,1);
   funct_Q2U(U,Q1);
   Boundaries(U);

   RK1D(U,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   Boundaries(U);
   
#elif DIM == 2 || DIM == 4

   Prim2Cons_All(Q,U);

   RK2D(U,Q,Q1,Q2,1);
   funct_Q2U(U,Q1);
   Boundaries(U);

   RK2D(U,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   Boundaries(U);
   
#elif DIM == 3 
   
   Prim2Cons_All(Q,U);
   RK3D(U,Q,Q1,Q2,1);
   funct_Q2U(U1,Q1);
   Boundaries(U1);
   RK3D(U1,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   Boundaries(U);
   
#endif

   return 0;
}
