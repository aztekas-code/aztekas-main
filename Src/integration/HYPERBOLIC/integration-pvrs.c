/**
 * @file integration-pvrs.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the primitive variables
 * \f$ \mathbf{U} \f$.
 */

//Do not erase any of these libraries//
#include"main.h"

void Integration()
{
   int n, i, j, k;

   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
#if DIM == 1 
   
   RK1D(U,U,U1,U2,1);
   Boundaries(U1);

   RK1D(U1,U,U1,U2,2);
   Boundaries(U2);
   U  = U2;
   U0 = U;
   
#elif DIM == 2 

   RK2D(U,U,U1,U2,1);
   Boundaries(U1);

   RK2D(U1,U,U1,U2,2);
   Boundaries(U2);
   U = U2;
   
#elif DIM == 3 
   
   funct_U2Q(Q,U);
   RK3D(U,Q,Q1,Q2,1);
   funct_Q2U(U1,Q1);
   Boundaries(U1);
   RK3D(U1,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   Boundaries(U);
   
#endif
}
