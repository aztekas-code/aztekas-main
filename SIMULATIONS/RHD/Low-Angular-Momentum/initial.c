#include"main.h"

void Initial()
{
   FILE *file;
   int idum;
   double dum;
   char line[100];
   double r;
   double t;
   double vr;
   double alpha, beta;
   double M = Black_Hole_Mass;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1 

   ////////////////////////////
   //-------Michel-1D--------//
   ////////////////////////////

   for(int i = 0; i <= Nx1; i++)
   {
      r = grid.X1[i];

      U(0,i) = density_0;//rho(r);
      U(1,i) = pressure_0;//pre(r);
      U(2,i) = 0.0;//velocity_0;
   }

#elif DIM == 2

   ////////////////////////////
   //-------Michel-2D--------//
   ////////////////////////////

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         r = grid.X1[i];

         U(0,i,j) = rho(r);
         U(1,i,j) = pre(r);
         U(2,i,j) = velocity_0;
         U(3,i,j) = 0.0;
      }
   }

#elif DIM == 4

   //////////////////////////////
   //-------Michel-2.5D--------//
   //////////////////////////////


   for(int j = 0; j <= Nx2; j++)
   {
      file = fopen("./Michel/analytic","r");
      for(int i = 0; i <= Nx1; i++)
      {

         r     = grid.X1[i];
         t     = grid.X2[j];
         alpha = sqrt(r/(r + 2.0*M));
         beta  = 2.0*M/r;

         idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,&U(RHO,i,j),&U(PRE,i,j),&U(VX1,i,j));
         vr    = (r/(r + 2.0*M))*U(VX1,i,j);

         U(VX2,i,j) = 0.0;
         if(r > 11.86)
         {
            U(VX3,i,j) = l_0*ftheta(t)*(alpha - beta*vr);
         }
         else
         {
            U(VX3,i,j) = 0.0;
         }
      }
      fclose(file);
   }

#endif
}
