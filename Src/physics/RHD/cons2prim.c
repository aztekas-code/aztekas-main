/**
 * @file /RHD/cons2prim.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Function that convert from Conservative to Primitives (RHD).
 *
 */

#include"main.h"
    
int Cons2Prim(double *u, double *q)
{
   int count;
   double D, tau, S_cov[3], S_con[3];
   double h, derh, f, derf;
   double Lorentz, SS;
   double theta, theta_0;
   gauge_ local_grid;
   
#if DIM == 1

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,tau,S_cov,S_con,h,derh,f,derf,\
                           Lorentz,SS,theta,theta_0,count,local_grid)
#endif
   for(int i = 0; i <= Nx1-0; i++)
   {
      local_grid.x[0] = grid.time;
      local_grid.x[1] = grid.X1[i];
      local_grid.x[2] = 0.0;
      local_grid.x[3] = 0.0;
      #if COORDINATES == SPHERICAL
      local_grid.x[2] = M_PI_2;
      #endif

      Get_Metric_Components(&local_grid);

      D        = q(DEN,i);
      tau      = q(ENE,i);
      S_cov[0] = q(MX1,i);
      S_cov[1] = 0.0;
      S_cov[2] = 0.0;

      S_con[0] = local_grid.gamma_con[0][0]*S_cov[0] + \
                 local_grid.gamma_con[0][1]*S_cov[1] + \
                 local_grid.gamma_con[0][2]*S_cov[2];
      S_con[1] = local_grid.gamma_con[1][0]*S_cov[0] + \
                 local_grid.gamma_con[1][1]*S_cov[1] + \
                 local_grid.gamma_con[1][2]*S_cov[2];
      S_con[2] = local_grid.gamma_con[2][0]*S_cov[0] + \
                 local_grid.gamma_con[2][1]*S_cov[1] + \
                 local_grid.gamma_con[2][2]*S_cov[2];
      
      SS = S_cov[0]*S_con[0] + S_cov[1]*S_con[1] + S_cov[2]*S_con[2];

      theta_0 = U(1,i)/U(0,i);
      f       = 1.0;
      count   = 0;

      while(fabs(f) > 0.00000001)
      {
      #if EOS == IDEAL
         h    = 1.0 + (K / (K - 1.0))*theta_0;
         derh = K / (K - 1.0);
      #elif EOS == DUST
         h    = 1.0;
         derh = 0.0;
      #elif EOS == STIFF
         h    = (K / (K - 1.0))*theta_0;
         derh = K / (K - 1.0);
      #endif

         Lorentz = sqrt(1.0 + SS/(D*D*h*h));

         f    = h*Lorentz - (theta_0/Lorentz) - (tau/D) - 1.0;
         derf = (1.0/Lorentz)*(derh - 1.0 - theta_0*((Lorentz*Lorentz - 1.0)/(Lorentz*Lorentz))*(derh/h));

         theta   = theta_0 - f/derf;
         theta_0 = theta;
         count++;

         if(count == 100000)
         {
            printf("                                          \n");
            printf("Spend too much time in Newton-Rhapson.\n");
            CHECK_NAN = TRUE;
            U = U0;
            Print_Values(&grid.time,&theta,&CHECK_NAN);
            exit(EXIT_FAILURE);
         }
      }

      #if EOS == IDEAL
      h    = 1.0 + (K / (K - 1.0))*theta_0;
      #elif EOS == DUST
      h    = 1.0;
      #elif EOS == STIFF
      h    = (K / (K - 1.0))*theta_0;
      #endif

      Lorentz = sqrt(1.0 + SS/(D*D*h*h));

      u(RHO,i) = D / Lorentz;
      u(PRE,i) = D*h*Lorentz - tau - D;
      u(VX1,i) = S_cov[0]/(D*h*Lorentz);
   }

#elif DIM == 2

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,tau,S_cov,S_con,h,derh,f,derf,\
                           Lorentz,SS,theta,theta_0,count,local_grid) \
                           collapse(2)
#endif
   for(int i = gc; i <= Nx1-gc; i++)
   {
      for(int j = gc; j <= Nx2-gc; j++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;
         #if POLAR == TRUE
         local_grid.x[2] = M_PI_2;
         #endif

         // Needs to be check. Fails in theta = Pi
         f = 1.0;
         if(x2max == M_PI && j == Nx2-gc){f = 0.0;}

         Get_Metric_Components(&local_grid);

         D        = q(DEN,i,j);
         tau      = q(ENE,i,j);
         S_cov[0] = q(MX1,i,j);
         S_cov[1] = q(MX2,i,j);
         S_cov[2] = 0.0;

         S_con[0] = local_grid.gamma_con[0][0]*S_cov[0] + \
                    local_grid.gamma_con[0][1]*S_cov[1] + \
                    local_grid.gamma_con[0][2]*S_cov[2];
         S_con[1] = local_grid.gamma_con[1][0]*S_cov[0] + \
                    local_grid.gamma_con[1][1]*S_cov[1] + \
                    local_grid.gamma_con[1][2]*S_cov[2];
         S_con[2] = local_grid.gamma_con[2][0]*S_cov[0] + \
                    local_grid.gamma_con[2][1]*S_cov[1] + \
                    local_grid.gamma_con[2][2]*S_cov[2];
         
         SS = S_cov[0]*S_con[0] + S_cov[1]*S_con[1] + S_cov[2]*S_con[2];

         theta_0 = U(1,i,j)/U(0,i,j);
         count   = 0;

         while(fabs(f) > 0.00000001)
         {
         #if EOS == IDEAL
            h    = 1.0 + (K / (K - 1.0))*theta_0;
            derh = K / (K - 1.0);
         #elif EOS == DUST
            h    = 1.0;
            derh = 0.0;
         #elif EOS == STIFF
            h    = (K / (K - 1.0))*theta_0;
            derh = K / (K - 1.0);
         #endif

            Lorentz = sqrt(1.0 + SS/(D*D*h*h));

            f    = h*Lorentz - (theta_0/Lorentz) - (tau/D) - 1.0;
            derf = (1.0/Lorentz)*(derh - 1.0 - theta_0*((Lorentz*Lorentz - 1.0)/(Lorentz*Lorentz))*(derh/h));

            theta   = theta_0 - f/derf;
            theta_0 = theta;
            count++;

            if(count == 100000)
            {
               printf("%e %e %e %e\n",Lorentz,h,derh,SS);
               printf("                                          \n");
               printf("Spend too much time in Newton-Rhapson.\n");
               CHECK_NAN = TRUE;
               U = U0;
               Print_Values(&grid.time,&theta,&CHECK_NAN);
               exit(EXIT_FAILURE);
            }
         }

         #if EOS == IDEAL
         h    = 1.0 + (K / (K - 1.0))*theta_0;
         #elif EOS == DUST
         h    = 1.0;
         #elif EOS == STIFF
         h    = (K / (K - 1.0))*theta_0;
         #endif

         Lorentz = sqrt(1.0 + SS/(D*D*h*h));

         u(RHO,i,j) = D / Lorentz;
         u(PRE,i,j) = D*h*Lorentz - tau - D;
         u(VX1,i,j) = S_cov[0]/(D*h*Lorentz);
         u(VX2,i,j) = S_cov[1]/(D*h*Lorentz);
      }
   }

#elif DIM == 4

#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp for private(D,tau,S_cov,S_con,h,derh,f,derf,\
                           Lorentz,SS,theta,theta_0,count,local_grid) \
                           collapse(2)
#endif
   for(int j = 0; j <= Nx2-0; j++)
   {
      for(int i = 0; i <= Nx1-0; i++)
      {
         local_grid.x[0] = grid.time;
         local_grid.x[1] = grid.X1[i];
         local_grid.x[2] = grid.X2[j];
         local_grid.x[3] = 0.0;

         Get_Metric_Components(&local_grid);

         D        = q(DEN,i,j);
         tau      = q(ENE,i,j);
         S_cov[0] = q(MX1,i,j);
         S_cov[1] = q(MX2,i,j);
         S_cov[2] = q(MX3,i,j);

         S_con[0] = local_grid.gamma_con[0][0]*S_cov[0] + \
                    local_grid.gamma_con[0][1]*S_cov[1] + \
                    local_grid.gamma_con[0][2]*S_cov[2];
         S_con[1] = local_grid.gamma_con[1][0]*S_cov[0] + \
                    local_grid.gamma_con[1][1]*S_cov[1] + \
                    local_grid.gamma_con[1][2]*S_cov[2];
         S_con[2] = local_grid.gamma_con[2][0]*S_cov[0] + \
                    local_grid.gamma_con[2][1]*S_cov[1] + \
                    local_grid.gamma_con[2][2]*S_cov[2];
         
         SS = S_cov[0]*S_con[0] + S_cov[1]*S_con[1] + S_cov[2]*S_con[2];

         theta_0 = U(PRE,i,j)/U(RHO,i,j);
         f       = 1.0;
         count   = 0;

         while(fabs(f) > 0.0000001)
         {
         #if EOS == IDEAL
            h    = 1.0 + (K / (K - 1.0))*theta_0;
            derh = K / (K - 1.0);
         #elif EOS == DUST
            h    = 1.0;
            derh = 0.0;
         #elif EOS == STIFF
            h    = (K / (K - 1.0))*theta_0;
            derh = K / (K - 1.0);
         #endif

            Lorentz = sqrt(1.0 + SS/(D*D*h*h));

            f    = h*Lorentz - (theta_0/Lorentz) - (tau/D) - 1.0;
            derf = (1.0/Lorentz)*(derh - 1.0 - theta_0*((Lorentz*Lorentz - 1.0)/(Lorentz*Lorentz))*(derh/h));

            theta   = theta_0 - f/derf;
            theta_0 = theta;
            count++;

            if(count == 10000000)
            {
               printf("                                          \n");
               printf("Spend too much time in Newton-Rhapson.\n");
               printf("%d %d %e %e %e %e\n",i,j,D,tau,SS,f);
               CHECK_NAN = TRUE;
               U = U0;
               Print_Values(&grid.time,&theta,&CHECK_NAN);
               exit(EXIT_FAILURE);
            }
         }

         #if EOS == IDEAL
         h    = 1.0 + (K / (K - 1.0))*theta_0;
         #elif EOS == DUST
         h    = 1.0;
         #elif EOS == STIFF
         h    = (K / (K - 1.0))*theta_0;
         #endif

         Lorentz = sqrt(1.0 + SS/(D*D*h*h));

         u(RHO,i,j) = D / Lorentz;
         u(PRE,i,j) = D*h*Lorentz - tau - D;
         u(VX1,i,j) = S_cov[0]/(D*h*Lorentz);
         u(VX2,i,j) = S_cov[1]/(D*h*Lorentz);
         u(VX3,i,j) = S_cov[2]/(D*h*Lorentz);
      }
   }

#elif DIM == 3

   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
      {
         for(int k = 0; k <= Nx3; k++)
         {
            D  = q(0,i,j,k);
            E  = q(1,i,j,k);
            S1 = q(2,i,j,k);
            S2 = q(3,i,j,k);
            S3 = q(4,i,j,k);
 
            u(0,i,j,k) = D;
            #if EOS == IDEAL
            u(1,i,j,k) = ((2.0*K-2.0)*D*E+(1.0-K)*pow(S3,2.0)+(1.0-K)*pow(S2,2.0)+(1.0-K)*pow(S1,2.0))/(2.0*D);
            #elif EOS == DUST
            u(1,i,j,k) = 0.0;
            #endif
            u(2,i,j,k) = S1/D;
            u(3,i,j,k) = S2/D;
            u(4,i,j,k) = S3/D;
         }
      }
   }

#endif
   
   return 0;
}
