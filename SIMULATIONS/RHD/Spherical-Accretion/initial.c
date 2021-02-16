/* 
 *  aztekas initial module
 *  Date of creation: 27-11-2019 10:44:40
 *  author: Alejandro Aguayo-Ortiz 
 */
#include"main.h"

void Initial()
{
   FILE *file;
   int idum;
   double dum;
   char line[100];
   double r;
   double t;
   double R, z;
   double c = 3e+10;

   //Initialize time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

#if DIM == 1

   file = fopen("./Michel/analytic.dat","r");
   for(int i = 0; i <= Nx1; i++)
   {
      r     = grid.X1[i];
      t     = M_PI_2;

      double a     = Black_Hole_Spin;
      double M     = Black_Hole_Mass;
      double rho2  = r*r + a*a*cos(t)*cos(t);
      double grr   = 1.0 + 2.0*M*r/rho2;
      double lapse = 1.0/sqrt(grr);
      double beta  = 2.0*M*r/rho2/(grr);

      U(RHO,i) = 1.0;
      //U(PRE,j) = ((K - 1.0)*pow(U(RHO,i,j),K)*cs*cs)/(K*(K - 1.0) - K*cs*cs)/pow(1.0/density_0,K-1.0);
      U(PRE,i) = Temp*U(RHO,i)*pow(U(RHO,i)/(1.0/density_0),K-1.0);
      U(VX1,i) = beta/lapse/grr;

      //idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,&U(RHO,i,j),&U(PRE,i,j),&U(VX1,i,j));
   }
   fclose(file);

#elif DIM == 2 || DIM == 4 

//   for(int j = 0; j <= Nx2; j++)
//   {
//      for(int i = 0; i <= Nx1; i++)
//      {
//         r     = grid.X1[i];
//         t     = grid.X2[j];
//
//         double a    = Black_Hole_Spin;
//         double M    = Black_Hole_Mass;
//         double rho2 = r*r + a*a*cos(t)*cos(t);
//         double gtt  = -(1.0 - 2.0*M*r/rho2);
//         double gtr  = 2.0*M*r/rho2;
//         double gtp  = -2.0*M*a*r*sin(t)*sin(t)/rho2;
//         double grr  = 1.0 + 2.0*M*r/rho2;
//         double grp  = -a*grr*sin(t)*sin(t);
//         double gpp  = sin(t)*sin(t)*(rho2 + a*a*grr*sin(t)*sin(t));
//         double gTT  = -grr;
//         double gTR  = gtr;
//
//         double den  = 1.0 - gTR*gtr;
//         double t1   = gTT/den;
//         double t2   = gTR*grr/den;
//         double t3   = gTR*grp/den;
//         double r1   = gtr*gTT/den;
//         double r2   = grr/den;
//         double r3   = grp/den;
//         double p1   = gtp*gTT;
//         double p2   = grp + gtp*gTR*grr/den;
//         double p3   = gpp + gtp*gTR*grp/den;
//
//         double UR;
//         double UP;
//         idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,&U(RHO,i,j),&U(PRE,i,j),&UR);
//         U(VX2,i,j) = 0.0;
//         UP = 0.0;
//
//         UP = UP + UR*a/(r*r - 2.0*M*r + a*a);
//
//         double a1 = t1;
//         double b1 = (t2 + r1)*UR + (t3 + p1)*UP;
//         double c1 = 1.0 + r2*UR*UR + (p2 + r3)*UR*UP + p3*UP*UP;
//
//         double Ut;
//         if(r < M + sqrt(M*M - a*a))
//         {
//            Ut = (-b1 + sqrt(b1*b1 - 4*a1*c1))/(2.0*a1);
//         }
//         else
//         {
//            Ut = (-b1 + sqrt(b1*b1 - 4*a1*c1))/(2.0*a1);
//         }
//
//         double Ur = r1*Ut + r2*UR + r3*UP;
//         double Up = p1*Ut + p2*UR + p3*UP;
//         double UT = t1*Ut + t2*UR + t3*UP;
//
//         double alpha = 1.0/sqrt(1.0 + 2.0*M*r/rho2);
//         double W     = alpha*UT;
//
//         U(VX1,i,j) = UR;///W;
//         U(VX3,i,j) = 0.0;
//      }
//      fclose(file);
//   }

//   density_0 = U(RHO,Nx1-gc,Nx2-gc);
   
   for(int j = 0; j <= Nx2; j++)
   {
      file = fopen("./Michel/analytic.dat","r");
      for(int i = 0; i <= Nx1; i++)
      {
         r     = grid.X1[i];
         t     = grid.X2[j];

         double a     = Black_Hole_Spin;
         double M     = Black_Hole_Mass;
         double rho2  = r*r + a*a*cos(t)*cos(t);
         double grr   = 1.0 + 2.0*M*r/rho2;
         double lapse = 1.0/sqrt(grr);
         double beta  = 2.0*M*r/rho2/(grr);

         U(RHO,i,j) = 1.0;
         //U(PRE,i,j) = ((K - 1.0)*pow(U(RHO,i,j),K)*cs*cs)/(K*(K - 1.0) - K*cs*cs)/pow(1.0/density_0,K-1.0);
         U(PRE,i,j) = Temp*U(RHO,i,j)*pow(density_0,K-1.0);
         U(VX1,i,j) = beta/lapse/grr;
         U(VX2,i,j) = 0.0;
         U(VX3,i,j) = 0.0;

//         idum = fscanf(file,"%lf %lf %lf %lf %lf\n",&dum,&U(RHO,i,j),&U(PRE,i,j),&U(VX1,i,j),&dum);
      }
      fclose(file);
   }
   
#endif
}
