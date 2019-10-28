#include"main.h"
    
void Prim2Cons_All(double *q, double *u)
{
   double E;
   double rho, p, vx1, vx2, vx3;
   double P[eq+1];
   eos_ eos;
   gauge_ local_grid;

#if DIM == 1

#ifdef _OPENMP
   #pragma omp parallel shared(grid)
   #pragma omp for private(rho,p,vx1,vx2,vx3,E,local_grid,eos,P)
#endif
   for(int i = 0; i <= Nx1; i++)
   {
      local_grid.x[0] = grid.time;
      local_grid.x[1] = grid.X1[i];
      local_grid.x[2] = 0.0;
      local_grid.x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      local_grid.x[2] = M_PI_2;
      #endif
    
      rho = u(0,i);
      p   = u(1,i);
      vx1 = u(2,i);
      vx2 = 0.0;
      vx3 = 0.0;

      P[0] = rho;
      P[1] = p;

      EoS(&eos,P,&local_grid);

      E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;

      q(RHO,i) = rho;
      q(PRE,i) = E;
      q(VX1,i) = rho*vx1;
   }

#elif DIM == 2

#if _OPEN_MP
   #pragma omp parallel shared (grid)
   #pragma omp for private(rho,p,vx1,vx2,vx3,E,local_grid,eos,P) collapse(2)
#endif
   for(int j = 0; j <= Nx2-0; j++)
   {
      for(int i = 0; i <= Nx1-0; i++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;
         #if POLAR == TRUE
         local_grid.x[2] = M_PI_2;
         #endif
       
         rho = u(0,i,j);
         p   = u(1,i,j);
         vx1 = u(2,i,j);
         vx2 = u(3,i,j);
         vx3 = 0.0;

         P[0] = rho;
         P[1] = p;

         EoS(&eos,P,&local_grid);

         E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
   
         q(RHO,i,j) = rho;
         q(PRE,i,j) = E;
         q(VX1,i,j) = rho*vx1;
         q(VX2,i,j) = rho*vx2;
      }
   }

#elif DIM == 4

#if _OPEN_MP
   #pragma omp parallel shared (grid)
   #pragma omp for private(rho,p,vx1,vx2,vx3,E,local_grid,eos,P) collapse(2)
#endif
   for(int i = 0; i <= Nx1-0; i++)
   {
      for(int j = 0; j <= Nx2-0; j++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;
       
         rho = u(0,i,j);
         p   = u(1,i,j);
         vx1 = u(2,i,j);
         vx2 = u(3,i,j);
         vx3 = u(4,i,j);

         P[0] = rho;
         P[1] = p;

         EoS(&eos,P,&local_grid);

         E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
   
         q(0,i,j) = rho;
         q(1,i,j) = E;
         q(2,i,j) = rho*vx1;
         q(3,i,j) = rho*vx2;
         q(4,i,j) = rho*vx3;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            local_grid.x[0] = grid.time;
            local_grid.x[1] = grid.X1[i];
            local_grid.x[2] = grid.X2[j];
            local_grid.x[3] = grid.X3[k];
        
            rho = u(0,i,j,k);
            p   = u(1,i,j,k);
            vx1 = u(2,i,j,k);
            vx2 = u(3,i,j,k);
            vx3 = u(4,i,j,k);

            P[0] = rho;
            P[1] = p;

            EoS(&eos,P,local_grid);

            E = 0.5 * rho * (vx1*vx1 + vx2*vx2 + vx3*vx3) + rho*eos.e;
          
            q(0,i,j,k) = rho;
            q(1,i,j,k) = E;
            q(2,i,j,k) = rho*vx1;
            q(3,i,j,k) = rho*vx2;
            q(4,i,j,k) = rho*vx3;
         }
      }
   }

#endif
}
