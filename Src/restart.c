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
#include"main.h"
#include"param.h"
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
