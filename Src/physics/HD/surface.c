#include"main.h"

void Surface_Volume()
{
   int i, j, k;

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      #if COORDINATES == CARTESIAN
      S1p(i) = 1.0;
      S1m(i) = 1.0;
      #elif COORDINATES == CYLINDRICAL
      S1p(i) = grid.X1p[i]/grid.X1[i];
      S1m(i) = grid.X1m[i]/grid.X1[i];
      #elif COORDINATES == SPHERICAL
      S1p(i) = 1.0;
      S1m(i) = 1.0;
      #endif
   }
   
#elif DIM == 2 || DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         #if COORDINATES == CARTESIAN
         S1p(i,j) = 1.0;
         S1m(i,j) = 1.0;
         S2p(i,j) = 1.0;
         S2m(i,j) = 1.0;
         #elif COORDINATES == CYLINDRICAL
         S1p(i,j) = grid.X1p[i]/grid.X1[i];
         S1m(i,j) = grid.X1m[i]/grid.X1[i];
         S2p(i,j) = 1.0;
         S2m(i,j) = 1.0;
         #elif COORDINATES == SPHERICAL
         S1p(i,j) = 1.0;
         S1m(i,j) = 1.0;
         S2p(i,j) = 1.0/grid.X1[i];
         S2m(i,j) = 1.0/grid.X1[i];
         #endif
         }
   }
#elif DIM == 3
   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            S1p(i,j,k) = 1.0;
            S1m(i,j,k) = 1.0;
            S2p(i,j,k) = 1.0;
            S2m(i,j,k) = 1.0;
            S3p(i,j,k) = 1.0;
            S3m(i,j,k) = 1.0;
         }
      }
   }

#endif
}
