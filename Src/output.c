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
#include<stdlib.h>
#include<string.h>
#include"./Headers/main.h"

int PrintValues(double *tprint, double *dtprint, int *itprint)
{
   int n, i, j, k;

   if(time >= *tprint)
   {
      printf("Time:%e, dt: %e, dx:%e \n", time, dt, dx1);

      if(graf == 1)
      {
         Output1(itprint);
      }
      else if(graf == 2)
      {
         Output2(itprint);
      }
      else if(graf == 3)
      {
         Output3(itprint);
      }

      *tprint = *tprint + *dtprint;
      ++*itprint;
   }

   time = time + dt;
   return 0;
}

int Output1(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"w");   

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"###################################\n");

   for(i = gc; i <= Nx1-gc; i++)
   {
      fprintf(file,"%e %e %e %e\n",X1[i],\
      U[c1(0,i)],\
      U[c1(1,i)],\
      U[c1(2,i)]);
   }

   printf("itprint %d \n",*itprint);
   fclose(file);

   return 0;
}

int Output2(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"w");

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
   fprintf(file,"###################################\n");

   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fprintf(file,"%e %e %e %e %e %e\n",X1[i],X2[j],U[c2(0,i,j)],\
         U[c2(1,i,j)],\
         U[c2(2,i,j)],\
         U[c2(3,i,j)]);
      }
   }

   printf("itprint %d \n",*itprint);
   fclose(file);

   return 0;
}

int Output3(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"w");

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
   fprintf(file,"%d \n",Nx3-2*gc+1);
   fprintf(file,"###################################\n");

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = Nx2/2; j <= Nx2/2; j++)
      {
         for(k = Nx3/2; k <= Nx3/2; k++)
         {
            fprintf(file,"%f %f %f %f %f %f %f %f\n",X1[i],X2[j],X3[k],\
            U[c3(0,i,j,k)],\
            U[c3(1,i,j,k)],\
            U[c3(2,i,j,k)],\
            U[c3(3,i,j,k)],\
            U[c3(4,i,j,k)]);
         }
      }
   }

   printf("itprint %d \n",*itprint);
   fclose(file);

   return 0;
}
