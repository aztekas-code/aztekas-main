//Do not erase any of these libraries//
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"./Headers/main.h"

//////////////////////////////////////
////////---Riemann-Problem---/////////
//////////////////////////////////////
void RIEMANN()
{
   int i, j, k, n; 

   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         if(X1[i] < x_0)
         {
            U[c1(0,i)] = nl;
            U[c1(1,i)] = pl;
            U[c1(2,i)] = vx1l;
         }
         else
         {
            U[c1(0,i)] = nr;
            U[c1(1,i)] = pr;
            U[c1(2,i)] = vx1r;
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(X1[i] < x_0)
            {
               U[c2(0,i,j)] = nl;
               U[c2(1,i,j)] = pl;
               U[c2(2,i,j)] = vx1l;
               U[c2(3,i,j)] = vx2l;
            }
            else
            {
               U[c2(0,i,j)] = nr;
               U[c2(1,i,j)] = pr;
               U[c2(2,i,j)] = vx1r;
               U[c2(3,i,j)] = vx2r;
            }
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            for(k = 0; k <= Nx3; k++)
            {
               if(X3[k] < x_0)
               {
                  U[c3(0,i,j,k)] = nl;
                  U[c3(1,i,j,k)] = pl;
                  U[c3(2,i,j,k)] = vx1l;
                  U[c3(3,i,j,k)] = vx2l;
                  U[c3(4,i,j,k)] = vx3l;
               } 
               else 
               { 
                  U[c3(0,i,j,k)] = nr;
                  U[c3(1,i,j,k)] = pr;
                  U[c3(2,i,j,k)] = vx1r;
                  U[c3(3,i,j,k)] = vx2r;
                  U[c3(4,i,j,k)] = vx3r;
               }
            }
         }
      }
   }
}

//////////////////////////////////////
////////---Kelvin-Helmholtz--/////////
//////////////////////////////////////
void KH()
{
   int i, j, k, n;

   if(dim == 1)
   {
      printf("Impossible to test Kelvin-Helmholtz's instability in 1D.\n");
      exit(1);
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(fabs(X2[j]) < x_0)
            {
               U[c2(0,i,j)] = nr;
               U[c2(1,i,j)] = pr;
               U[c2(2,i,j)] = vx1r*(1 + 0.01*cos(1*M_PI*X1[i])*cos(1*M_PI*X2[j]));
               U[c2(3,i,j)] = vx2r*(1 + 0.01*cos(1*M_PI*X1[i])*cos(1*M_PI*X2[j]));
            }
            else
            {
               U[c2(0,i,j)] = nl;
               U[c2(1,i,j)] = pl;
               U[c2(2,i,j)] = vx1l*(1 + 0.01*cos(1*M_PI*X1[i])*cos(1*M_PI*X2[j]));
               U[c2(3,i,j)] = vx2l*(1 + 0.01*cos(1*M_PI*X1[i])*cos(1*M_PI*X2[j]));
            }
         }
      }
   }
}

//////////////////////////////////////
////////---------Jet---------/////////
//////////////////////////////////////
void JET()
{
   int i, j, k, n;
   
   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         if(fabs(X1[i]) <= r_jet)
         {
            U[c1(0,i)] = n_jet;
            U[c1(1,i)] = p_jet;
            U[c1(2,i)] = vx1_jet;
         }
         else
         {
            U[c1(0,i)] = n_atm;
            U[c1(1,i)] = p_atm;
            U[c1(2,i)] = vx1_atm;
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(fabs(X1[i]) <= r_jet)
            {
               if(fabs(X2[j]) <= z_jet)
               {
                  U[c2(0,i,j)] = n_jet;
                  U[c2(1,i,j)] = p_jet;
                  U[c2(2,i,j)] = vx1_jet;
                  U[c2(3,i,j)] = vx2_jet;
               }
               else
               {
                  U[c2(0,i,j)] = n_atm;
                  U[c2(1,i,j)] = p_atm;
                  U[c2(2,i,j)] = vx1_atm;
                  U[c2(3,i,j)] = vx2_atm;
               }
            }
            else
            {
               U[c2(0,i,j)] = n_atm;
               U[c2(1,i,j)] = p_atm;
               U[c2(2,i,j)] = vx1_atm;
               U[c2(3,i,j)] = vx2_atm;
            }
         }
      }
   }
}

//////////////////////////////////////
////////-Spherical-Accretion-/////////
//////////////////////////////////////
void SPH_ACC()
{
   int i, j, k, n;
   double x1, x2, x3, r;

   if(dim == 1)
   {
      for(i = 0; i <= Nx1; i++)
      {
         r = X1[i];

         U[c1(0,i)] = density_0;
         U[c1(1,i)] = pressure_0;
         U[c1(2,i)] = velocity_0;
      
         if(r > r_out)
         {
            U[c1(2,i)] = velocity_0;
         }
      }
   }
   else if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)
      {   
         for(j = 0; j <= Nx2; j++)
         {
            if (alfa <= 1)
            {
               x1 = X1[i];
               x2 = X2[j];
               r = sqrt(x1*x1 + x2*x2);

               U[c2(0,i,j)] = density_0;
               U[c2(1,i,j)] = pressure_0;
               U[c2(2,i,j)] = velocity_0*(x1/r);
               U[c2(3,i,j)] = velocity_0*(x2/r);

               if(r > r_out)
               {
                  U[c2(2,i,j)] = velocity_0*(x1/r);
                  U[c2(3,i,j)] = velocity_0*(x2/r);
               }
            }
            else if (alfa == 2)
            {
               x1 = X1[i];
               x2 = X2[j];
               r = x1;

               U[c2(0,i,j)] = density_0;
               U[c2(1,i,j)] = pressure_0;
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = 0.0;

               if(r > r_out)
               {
                  U[c2(2,i,j)] = velocity_0*(x1/r);
               }
            }
         }
      }
   }
}

//////////////////////////////////////
////////---------Wind--------/////////
//////////////////////////////////////
void WIND()
{
   int i, j, k;
   double r;

   if(dim == 2)
   {
      for(i = 0; i <= Nx1; i++)                                                 
      {                                                                         
         for(j = 0; j <= Nx2; j++)                                              
         {                                                                      
            if(alfa <= 1)
            {
               U[c2(0,i,j)] = density_0;                                              
               U[c2(1,i,j)] = pressure_0;
               U[c2(2,i,j)] = 0.0;                                              
               U[c2(3,i,j)] = velocity_0;
            }
            else if(alfa == 2)
            {
               U[c2(0,i,j)] = density_0;                                     
               U[c2(1,i,j)] = pressure_0;                                    
               U[c2(2,i,j)] =  velocity_0*cos(X2[j])*sqrt(1 + 2/X1[i]);                                           
               U[c2(3,i,j)] = -velocity_0*sin(X2[j])*X1[i];                                    
            }
          }                                                                      
       }                                          
   }
}
