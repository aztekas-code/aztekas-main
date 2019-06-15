/**
 * @file bound_cond.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Standard boundary conditions.
 */

//Do not erase any of these libraries//
#include"main.h"

void Outflow(double *B)
{
   int i, j, k, n, cell;
   
#if DIM == 1

   //////////////////////////////////////////////////////////
   ///////////////////// Outflow-1D /////////////////////////
   //////////////////////////////////////////////////////////

   for(cell = 0; cell < gc; cell++)
   {
      for(n = 0; n < eq; n++)
      {
         // r = 0 boundary in cylindrical and spherical coordinates
      #if alfa > 0 && outflow_x1min == 1
         if(x1min == 0.0)
            B[c1(n,gc)] = B[c1(n,gc+1)];
      #endif

      #if outflow_x1max == 1
         B[c1(n,Nx1-cell)] = B[c1(n,Nx1-gc)]; //x1max
      #endif
      #if outflow_x1min == 1
         B[c1(n,cell)] = B[c1(n,gc)]; //x1min
      #endif
      }
   }
   //////////////////////////////////////////////////////////

#elif DIM == 2 || DIM == 4

   //////////////////////////////////////////////////////////
   ///////////////////// Outflow21D /////////////////////////
   //////////////////////////////////////////////////////////

   // Outflow in X1 //
   for(j = 0; j <= Nx2; j++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(n = 0; n < eq; n++)
         {
         // r = 0 boundary in cylindrical and spherical coordinates
         #if alfa >= 1 && outflow_x1min == 1
            if(x1min == 0.0) 
                B[c2(n,gc,j)] = B[c2(n,gc+1,j)];
         #endif

         #if outflow_x1max == 1
            B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-gc,j)]; //x1max
         #endif
         #if outflow_x1min == 1
            B[c2(n,cell,j)] = B[c2(n,gc,j)]; //x1min
         #endif
         }
      }
   }
   ///////////////////

   ///////////////////
   // Outflow in X2 //
   ///////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(n = 0; n < eq; n++)
         {
         // theta = 0 or M_PI boundary in spherical coordinates
         #if alfa == 2 && outflow_x2max == 1
            if(x2max/M_PI == 1.0)
               B[c2(n,i,Nx2-gc)] = B[c2(n,i,Nx2-gc-1)]; //x2max
         #endif
         #if alfa == 2 && outflow_x2max == 1
            if(x2min/M_PI == 0.0)
               B[c2(n,i,gc)] = B[c2(n,i,gc+1)]; //x2min
         #endif

         #if outflow_x2max == 1
            B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-gc)]; //x2max
         #endif
         #if outflow_x2min == 1
            B[c2(n,i,cell)] = B[c2(n,i,gc)]; //x2min
         #endif
         }
      }
   }
   //////////////////////////////////////////////////////////

#endif
}

void Reflection(double *B)
{
   int i, j, k, n, cell;

#if DIM == 1

   ///////////////////////////////////////////////////////////
   /////////////////////// Reflection-1D /////////////////////
   ///////////////////////////////////////////////////////////
   
   //////////////////////
   // Reflection on X1 //
   //////////////////////
   for(cell == 0; cell < gc; cell++)
   {
   // Reflection on x1max //
   #if reflective_x1max == 1
      for(n = 0; n < eq; n++)
      {
         B[c1(n,Nx1-cell)] = B[c1(n,Nx1-2*gc+cell+1)];
      }

      B[c1(2,Nx1-cell)] = -B[c1(2,Nx1-2*gc+cell+1)];
      
   #endif

   // Reflection on x1min //
   #if reflective_x1min == 1
      for(n = 0; n < eq; n++)
      {
      #if alfa == 0 
         B[c1(n,cell)] = B[c1(n,2*gc-cell-1)];
      #elif alfa >= 1
         B[c1(n,cell)] = B[c1(n,2*gc-cell)];
         B[c1(n,gc)] = B[c1(n,gc+1)];
      #endif
      }

      #if alfa == 0 
      B[c1(2,cell)] = -B[c1(2,2*gc-cell-1)];
      #elif alfa >= 1
      B[c1(2,cell)] = -B[c1(2,2*gc-cell)];
      B[c1(2,gc)] = 0.0;
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
   for(j = 0; j <= Nx2; j++)
   {
      for(cell = 0; cell < gc; cell++)
      {
      // Reflection on x1max //
      #if reflective_x1max == 1
         for(n = 0; n < eq; n++)
         {
            B[c2(n,Nx1-cell,j)] = B[c2(n,Nx1-2*gc+cell+1,j)];
         }

         B[c2(2,Nx1-cell,j)] = -B[c2(2,Nx1-2*gc+cell+1,j)];
      #endif

      // Reflection on x1min //
      #if reflective_x1min == 1
         for(n = 0; n < eq; n++)
         {
         #if alfa == 0 
            B[c2(n,cell,j)] = B[c2(n,2*gc-cell-1,j)];
         #elif alfa >= 1
            B[c2(n,cell,j)] = B[c2(n,2*gc-cell,j)];
            B[c2(n,gc,j)] = B[c2(n,gc+1,j)];
         #endif
         }

         #if alfa == 0 
         B[c2(2,cell,j)] = -B[c2(2,2*gc-cell-1,j)];
         #elif alfa >= 1
         B[c2(2,cell,j)] = -B[c2(2,2*gc-cell,j)];
         B[c2(2,gc,j)] = 0.0;
         #endif

      #endif
      }
   }
   //////////////////////

   //////////////////////
   // Reflective on X2 //
   //////////////////////
   for(i = 0; i <= Nx1; i++)
   {
      for(cell = 0; cell < gc; cell++)
      {
      // Reflection on x2max //
      #if reflective_x2max == 1
         for(n = 0; n < eq; n++)
         {
         #if alfa <= 1 
            B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-2*gc+cell+1)];
         #elif alfa == 2
            if(fabs(x2max - M_PI) <= 1.0e-05 || fabs(x2max - M_PI_2) <= 1.0e-05)
            {
               B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-2*gc+cell)];
               B[c2(n,i,Nx2-gc)] = B[c2(n,i,Nx2-gc-1)];
            }
            else
            {
               B[c2(n,i,Nx2-cell)] = B[c2(n,i,Nx2-2*gc+cell+1)];
             
            }
         #endif
         }

         #if alfa <= 1 
         B[c2(3,i,Nx2-cell)] = -B[c2(3,i,Nx2-2*gc+cell+1)];
         #elif alfa == 2
         if(fabs(x2max - M_PI) <= 1.0e-05 || fabs(x2max - M_PI_2) <= 1.0e-05)
         {
            B[c2(3,i,Nx2-cell)] = -B[c2(3,i,Nx2-2*gc+cell)];
            B[c2(3,i,Nx2-gc)] = 0.0;
         }
         else
         {
            B[c2(3,i,Nx2-cell)] = -B[c2(3,i,Nx2-2*gc+cell+1)];
         }
         #endif
      #endif

      // Reflection on x2min //
      #if reflective_x2min == 1
         for(n = 0; n < eq; n++)
         {
         #if alfa <= 1 
            B[c2(n,i,cell)] = B[c2(n,i,2*gc-cell-1)];
         #elif alfa == 2 
            if(fabs(x2min) <= 1.0e-05)
            {
               B[c2(n,i,cell)] = B[c2(n,j,2*gc-cell)];
               B[c2(n,i,gc)] = B[c2(n,i,gc+1)];
            }
            else
            {
               B[c2(n,i,cell)] = B[c2(n,i,2*gc-cell-1)];
            }
         #endif
         }

         #if alfa <= 1 
         B[c2(3,i,cell)] = -B[c2(3,i,2*gc-cell-1)];
         #elif alfa == 2
         if(fabs(x2min) <= 1.0e-05)
         {
            B[c2(3,i,cell)] = -B[c2(3,i,2*gc-cell)];
            B[c2(3,i,gc)] = 0.0;
         }
         else
         {
            B[c2(3,i,cell)] = -B[c2(3,i,2*gc-cell-1)];
         }
         #endif
      #endif
      }
   }
   //////////////////////
//////////////////////////////////////////////////////////////

#endif
}

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
  
   #if periodic_x1 == 1
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         B[c1(n,Nx1-cell)] = B[c1(n,2*gc-cell)];
         B[c1(n,cell)] = B[c1(n,Nx1-2*gc+cell)];
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
  
   #if periodic_x1 == 1
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            B[c2(n,Nx1-cell,j)] = B[c2(n,2*gc-cell,j)];
            B[c2(n,cell,j)] = B[c2(n,Nx1-2*gc+cell,j)];
         }
      }
   }
   #endif
   ////////////////////

   ////////////////////
   // Periodic on X2 //
   ////////////////////

   #if periodic_x2 == 1
   for(n = 0; n < eq; n++)
   {
      for(cell = 0; cell <= gc; cell++)
      {
         for(i = 0; i <= Nx1; i++)
         {
            B[c2(n,i,Nx2-cell)] = B[c2(n,i,2*gc-cell)];
            B[c2(n,i,cell)] = B[c2(n,i,Nx2-2*gc+cell)];
         }
      }
   }
   #endif
   ////////////////////
   //////////////////////////////////////////////////////////

#endif
}
