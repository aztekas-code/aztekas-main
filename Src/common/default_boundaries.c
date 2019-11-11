/**
 * @file default_boundaries.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Standard boundary conditions. Outflow, Periodic and Reflection.
 */

#include"main.h"

/*!
 * The function \b Outflow(), receives the vector solution
 * as an parameter \b B. It fills the value of the ghost cells in the specified
 * direction using the value of the last computed cell of the domain.
 */
void Outflow(double *B)
{
#if DIM == 1

   //////////////////////////////////////////////////////////
   ///////////////////// Outflow-1D /////////////////////////
   //////////////////////////////////////////////////////////

   #pragma omp parallel
   {
      #pragma omp for collapse(2)
      for(int cell = 0; cell < gc; cell++)
      {
         for(int n = 0; n < eq; n++)
         {
            // r = 0 boundary in cylindrical and spherical coordinates
         #if COORDINATES != CARTESIAN && outflow_x1min == TRUE
            if(x1min == 0.0)
               B(n,gc) = B(n,gc+1);
         #endif

         #if outflow_x1max == TRUE
            B(n,Nx1-cell) = B(n,Nx1-gc); //x1max
         #endif
         #if outflow_x1min == TRUE
            B(n,cell) = B(n,gc); //x1min
         #endif
         }
      }
   }
   //////////////////////////////////////////////////////////

#elif DIM == 2 || DIM == 4

   //////////////////////////////////////////////////////////
   ///////////////////// Outflow21D /////////////////////////
   //////////////////////////////////////////////////////////

   // Outflow in X1 //
   #pragma omp parallel
   {
      #pragma omp for collapse(3)
      for(int j = 0; j <= Nx2; j++)
      {
         for(int cell = 0; cell <= gc; cell++)
         {
            for(int n = 0; n < eq; n++)
            {
            // r = 0 boundary in cylindrical and spherical coordinates
            #if COORDINATES != CARTESIAN && outflow_x1min == TRUE
               if(x1min == 0.0) 
                   B(n,gc,j) = B(n,gc+1,j);
            #endif

            #if outflow_x1max == TRUE
               B(n,Nx1-cell,j) = B(n,Nx1-gc,j); //x1max
            #endif
            #if outflow_x1min == TRUE
               B(n,cell,j) = B(n,gc,j); //x1min
            #endif
            }
         }
      }
   }
   ///////////////////

   ///////////////////
   // Outflow in X2 //
   ///////////////////
   #pragma omp parallel
   {
      #pragma omp for collapse(3)
      for(int cell = 0; cell <= gc; cell++)
      {
         for(int i = 0; i <= Nx1; i++)
         {
            for(int n = 0; n < eq; n++)
            {
            // theta = 0 or M_PI boundary in spherical coordinates
            #if COORDINATES == SPHERICAL && outflow_x2max == TRUE
               if(x2max/M_PI == 1.0)
                  B(n,i,Nx2-gc) = B(n,i,Nx2-gc-1); //x2max
            #endif
            #if COORDINATES == SPHERICAL && outflow_x2max == TRUE
               if(x2min/M_PI == 0.0)
                  B(n,i,gc) = B(n,i,gc+1); //x2min
            #endif

            #if outflow_x2max == TRUE
               B(n,i,Nx2-cell) = B(n,i,Nx2-gc); //x2max
            #endif
            #if outflow_x2min == TRUE
               B(n,i,cell) = B(n,i,gc); //x2min
            #endif
            }
         }
      }
   }
   //////////////////////////////////////////////////////////

#endif
}

/*!
 * The function \b Reflection(), receives the vector solution
 * as an parameter \b B. It fills the value of the ghost cells in the specified
 * direction using the value of the mirrored cells, and for the velocity it
 * changes sign.
 */
void Reflection(double *B)
{
#if DIM == 1

   ///////////////////////////////////////////////////////////
   /////////////////////// Reflection-1D /////////////////////
   ///////////////////////////////////////////////////////////
   
   //////////////////////
   // Reflection on X1 //
   //////////////////////
   for(int cell = 0; cell < gc; cell++)
   {
   // Reflection on x1max //
   #if reflective_x1max == TRUE
      for(int n = 0; n < eq; n++)
      {
         B(n,Nx1-cell) = B(n,Nx1-2*gc+cell+1);
      }

      B(VX1,Nx1-cell) = -B(VX1,Nx1-2*gc+cell+1);
      
   #endif

   // Reflection on x1min //
   #if reflective_x1min == TRUE
      for(int n = 0; n < eq; n++)
      {
      #if COORDINATES == CARTESIAN
         B(n,cell) = B(n,2*gc-cell-1);
      #elif COORDINATES != CARTESIAN
         B(n,cell) = B(n,2*gc-cell);
         B(n,gc) = B(n,gc+1);
      #endif
      }

      #if COORDINATES == CARTESIAN
      B(VX1,cell) = -B(VX1,2*gc-cell-1);
      #elif COORDINATES != CARTESIAN
      B(VX1,cell) = -B(VX1,2*gc-cell);
      B(VX1,gc) = 0.0;
      #endif

   #endif
   }
   //////////////////////
   ///////////////////////////////////////////////////////////
 
#elif DIM == 2 || DIM == 4

   ///////////////////////////////////////////////////////////
   /////////////////////// Reflection-2D /////////////////////
   ///////////////////////////////////////////////////////////

   //////////////////////
   // Reflection on X1 //
   //////////////////////
   for(int j = 0; j <= Nx2; j++)
   {
      for(int cell = 0; cell < gc; cell++)
      {
      // Reflection on x1max //
      #if reflective_x1max == TRUE
         for(int n = 0; n < eq; n++)
         {
            B(n,Nx1-cell,j) = B(n,Nx1-2*gc+cell+1,j);
         }

         B(VX1,Nx1-cell,j) = -B(VX1,Nx1-2*gc+cell+1,j);
      #endif

      // Reflection on x1min //
      #if reflective_x1min == TRUE
         for(int n = 0; n < eq; n++)
         {
         #if COORDINATES == CARTESIAN
            B(n,cell,j) = B(n,2*gc-cell-1,j);
         #elif COORDINATES != CARTESIAN
            B(n,cell,j) = B(n,2*gc-cell,j);
            B(n,gc,j) = B(n,gc+1,j);
         #endif
         }

         #if COORDINATES == CARTESIAN
         B(VX1,cell,j) = -B(VX1,2*gc-cell-1,j);
         #elif COORDINATES != CARTESIAN
         B(VX1,cell,j) = -B(VX1,2*gc-cell,j);
         B(VX1,gc,j) = 0.0;
         #endif

      #endif
      }
   }
   //////////////////////

   //////////////////////
   // Reflective on X2 //
   //////////////////////
   for(int i = gc; i <= Nx1-gc; i++)
   {
      for(int cell = 0; cell < gc; cell++)
      {
      // Reflection on x2max //
      #if reflective_x2max == TRUE
         for(int n = 0; n < eq; n++)
         {
         #if COORDINATES != SPHERICAL
            B(n,i,Nx2-cell) = B(n,i,Nx2-2*gc+cell+1);
         #elif COORDINATES == SPHERICAL
            if(fabs(x2max - M_PI) <= 1.0e-05 || fabs(x2max - M_PI_2) <= 1.0e-05)
            {
               B(n,i,Nx2-cell) = B(n,i,Nx2-2*gc+cell);
               B(n,i,Nx2-gc) = B(n,i,Nx2-gc-1);
            }
            else
            {
               B(n,i,Nx2-cell) = B(n,i,Nx2-2*gc+cell+1);
             
            }
         #endif
         }

         #if COORDINATES != SPHERICAL
         B(VX2,i,Nx2-cell) = -B(VX2,i,Nx2-2*gc+cell+1);
         #elif COORDINATES == SPHERICAL
         if(fabs(x2max - M_PI) <= 1.0e-05 || fabs(x2max - M_PI_2) <= 1.0e-05)
         {
            B(VX2,i,Nx2-cell) = -B(VX2,i,Nx2-2*gc+cell);
            B(VX2,i,Nx2-gc)   =  0.0;
         }
         else
         {
            B(VX2,i,Nx2-cell) = -B(VX2,i,Nx2-2*gc+cell+1);
         }
         #endif
      #endif

      // Reflection on x2min //
      #if reflective_x2min == TRUE
         for(int n = 0; n < eq; n++)
         {
         #if COORDINATES != SPHERICAL
            B(n,i,cell) = B(n,i,2*gc-cell-1);
         #elif COORDINATES == SPHERICAL
            if(fabs(x2min) <= 1.0e-05)
            {
               B(n,i,cell) = B(n,i,2*gc-cell);
               B(n,i,gc)   = B(n,i,gc+1);
            }
            else
            {
               B(n,i,cell) = B(n,i,2*gc-cell-1);
            }
         #endif
         }

         #if COORDINATES != SPHERICAL
         B(VX2,i,cell) = -B(VX2,i,2*gc-cell-1);
         #elif COORDINATES == SPHERICAL
         if(fabs(x2min) <= 1.0e-05)
         {
            B(VX2,i,cell) = -B(VX2,i,2*gc-cell);
            B(VX2,i,gc)   =  0.0;
         }
         else
         {
            B(VX2,i,cell) = -B(VX2,i,2*gc-cell-1);
         }
         #endif
      #endif
      }
   }
   //////////////////////
//////////////////////////////////////////////////////////////

#endif
}

/*!
 * The function \b Periodic(), receives the vector solution
 * as an parameter \b B. It fills the value of the ghost cells with the values
 * of the correspondent other side of the domain.
 */
void Periodic(double *B)
{
   int i, j, k, n, cell;

#if DIM == 1

   //////////////////////////////////////////////////////////
   //////////////////// Periodic-1D /////////////////////////
   //////////////////////////////////////////////////////////

   ////////////////////
   // Periodic on X1 //
   ////////////////////
  
   #if periodic_x1 == TRUE
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         B(n,Nx1-cell) = B(n,2*gc-cell);
         B(n,cell) = B(n,Nx1-2*gc+cell);
      }
   }
   #endif
   //////////////////////////////////////////////////////////

#elif DIM == 2 || DIM == 4

   //////////////////////////////////////////////////////////
   //////////////////// Periodic-2D /////////////////////////
   //////////////////////////////////////////////////////////

   ////////////////////
   // Periodic on X1 //
   ////////////////////
  
   #if periodic_x1 == TRUE
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            B(n,Nx1-cell,j) = B(n,2*gc-cell,j);
            B(n,cell,j) = B(n,Nx1-2*gc+cell,j);
         }
      }
   }
   #endif
   ////////////////////

   ////////////////////
   // Periodic on X2 //
   ////////////////////

   #if periodic_x2 == TRUE
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(i = 0; i <= Nx1; i++)
         {
            B(n,i,Nx2-cell) = B(n,i,2*gc-cell);
            B(n,i,cell) = B(n,i,Nx2-2*gc+cell);
         }
      }
   }
   #endif
   ////////////////////
   //////////////////////////////////////////////////////////

#endif
}
