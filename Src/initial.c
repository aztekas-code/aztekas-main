/*
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<string.h>
#include"./Headers/main.h"

void INITIAL(double *dtprint)
{
   int n, i, j, k;

   //Initialize time
   time = 0.0;

   //Initialize dt
   dt = 0.0;

   switch(TEST)
   {
      case 0:
         printf("Custom problem in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         INIT_CUSTOM();
      break;

      case 1: //Riemann problem
         printf("Testing Riemann problem in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         RIEMANN();
      break;

      case 2: //Kelvin-Helmholtz instability
         printf("Testing Kelvin-Helmholtz instability problem in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         KH();
      break;

      case 3: //Jet
         printf("Testing Jet in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         JET();
      break;

      case 4:
         printf("Testing spherical accretion in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         SPH_ACC();
      break;

      case 5:
         printf("Testing wind in %dD\n",dim);
         printf("Press enter to continue\n");
         getchar();
         WIND();
      break;
   }

      ///////////////////////////
      //-----BHL-rel-----------//
      ///////////////////////////
/*     
      for(i = 0; i <= Nx1; i++)
      {   
         for(j = 0; j <= Nx2; j++)
         {
            double RR = sqrt(X1[i]*X1[i] + X2[j]*X2[j]);
            double R  = X1[i];
            double z  = X2[i];
            double v  = 0.5;
            double cs = 0.1;
            double r  = MM/(cs*cs + v*v);
            double rr = 1.1*MM;

            if(RR <= rr)
            {
               U[c2(0,i,j)] = 0.1;
               U[c2(1,i,j)] = cs*cs*U[c2(0,i,j)]*(K-1)/(K*(K-1) - cs*cs*K);
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = v;
            }
            else
            {
               U[c2(0,i,j)] = 1.0;
               U[c2(1,i,j)] = cs*cs*U[c2(0,i,j)]*(K-1)/(K*(K-1) - cs*cs*K);
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = v;
            }
         }
      }
*/
      /////////////////////////////


      ///////////////////////////
      //--------Sedov----------//
      ///////////////////////////
      /*
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(pow(X1[i] - 0.5,2.0) + pow(X2[j] - 0.5,2.0) < pow(0.15,2.0))
            {
               U[c2(0,i,j)] = 1.0;
               U[c2(1,i,j)] = 1.0;
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = 0.0;
            }
            else
            {
               U[c2(0,i,j)] = 1.0;
               U[c2(1,i,j)] = 0.1;
               U[c2(2,i,j)] = 0.0;
               U[c2(3,i,j)] = 0.0;
            }
         }
      }
      */
      /////////////////////////////

      ///////////////////////////
      //-------Riemann---------//
      ///////////////////////////
/*
      for(i = 0; i <= Nx1; i++)
      {
         for(j = 0; j <= Nx2; j++)
         {
            if(X2[j] < x_0)
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
*/
   /////////////////////////////
}

void RESTART()
{
   FILE *file;
   int i, j, k, idum;
   double dum;
   char line[100];

   //Initialize dt
   dt = 0.1;

   file = fopen(restartfile,"r");

   // Skip first line
   idum = fscanf(file,"%s\n",line) ;

   // Read time
   idum = fscanf(file,"%lf\n",&time) ;

   // Read Nx1
   idum = fscanf(file,"%d\n",&dum) ;

   // Read Nx2
   idum = fscanf(file,"%d\n",&dum) ;

   printf("%i %f\n",idum,time) ;

   // Skip third line   
   idum = fscanf(file,"%s\n",line) ;

   // Read rest of file an initialize variables      
   if(dim == 1)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,\
         &U[c1(0,i)],&U[c1(1,i)],&U[c1(2,i)]);
      }
   }
   else if(dim == 2)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            idum = fscanf(file,"%lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
            &U[c2(0,i,j)],&U[c2(1,i,j)],&U[c2(2,i,j)],&U[c2(3,i,j)]);
         }
      }
   }
   else if(dim == 3)
   {
      for(i = gc; i <= Nx1-gc; i++)
      {
         for(j = gc; j <= Nx2-gc; j++)
         {
            for(k = gc; k <= Nx3-gc; k++)
            {
               idum = fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf\n",\
               &dum,&dum,&dum,\
               &U[c3(0,i,j,k)],&U[c3(1,i,j,k)],\
               &U[c3(2,i,j,k)],&U[c3(3,i,j,k)],&U[c3(4,i,j,k)]);
            }
         }
      }
    }

   // CALLING BOUNDARIES TO GET GHOST CELLS RIGHT
    BOUNDARIES(U);

    fclose(file);
}                                        

