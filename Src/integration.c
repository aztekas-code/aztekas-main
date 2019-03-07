/**
 * @file integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"
#include"vector.h"

int INTEGRATION()
{
   int n, i, j, k;
   double uu[eq+1];
   double qq[eq+1];

   //Runge-Kutta 2th-Order and Piecewie Polynomial Reconstruction
#if dim == 1 
   
   funct_U2Q(Q,U);

   RK1D(U,Q,Q1,Q2,1);
   funct_Q2U(U,Q1);
   BOUNDARIES(U);

   RK1D(U,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   BOUNDARIES(U);
   
#elif dim == 2 || dim == 4

   funct_U2Q(Q,U);

   RK2D(U,Q,Q1,Q2,1);
   funct_Q2U(U,Q1);
   BOUNDARIES(U);

   RK2D(U,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   BOUNDARIES(U);
   
#elif dim == 3 
   
   funct_U2Q(Q,U);
   RK3D(U,Q,Q1,Q2,1);
   funct_Q2U(U1,Q1);
   BOUNDARIES(U1);
   RK3D(U1,Q,Q1,Q2,2);
   funct_Q2U(U,Q2);
   BOUNDARIES(U);
   
#endif

   return 0;
}
