/**
 * @file output.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Output functions: ASCII and Binary.
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"main.h"

int PrintValues(double *tprint, double *dtprint, int *itprint)
{
   int n, i, j, k;

   if(time >= *tprint)
   {
      printf("Time:%e, dt: %e, dx:%e \n", time, dt, dx1);

      if(graf == 1)
      {
         if(binary == 1)
         {
            Output1_bin(itprint);
         }
         else
         {
            Output1(itprint);
         }
      }
      else if(graf == 2)
      {
         if(binary == 1)
         {
            Output2_bin(itprint);
         }
         else
         {
            Output2(itprint);
         }
      }
      else if(graf == 3)
      {
         if(binary == 1)
         {
            Output3_bin(itprint);
         }
         else
         {
            Output3(itprint);
         }
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

int Output1_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;
   
   size_X1 = Nx1-2*gc+1;
   
   num = *itprint;
   strcpy(ext,".bin");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"wb");

   fwrite(&time, sizeof time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&X1[gc], sizeof X1[gc], 1, file);
   fwrite(&X1[Nx1-gc], sizeof X1[gc], 1, file);   
         
   for(i = gc; i <= Nx1-gc; i++)
	{
      fwrite(&U[c1(0,i)], sizeof U[c1(0,i)], 1, file);
      fwrite(&U[c1(1,i)], sizeof U[c1(1,i)], 1, file);
      fwrite(&U[c1(2,i)], sizeof U[c1(2,i)], 1, file);
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

int Output2_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;
   
   size_X1 = Nx1-2*gc+1;
   size_X2 = Nx2-2*gc+1;
   
   num = *itprint;
   strcpy(ext,".bin");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"wb");

   fwrite(&time, sizeof time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&X1[gc], sizeof X1[gc], 1, file);
   fwrite(&X1[Nx1-gc], sizeof X1[gc], 1, file);   
   fwrite(&X2[gc], sizeof X2[gc], 1, file);
   fwrite(&X2[Nx2-gc], sizeof X2[gc], 1, file);  
         
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fwrite(&U[c2(0,i,j)], sizeof U[c2(0,i,j)], 1, file);
         fwrite(&U[c2(1,i,j)], sizeof U[c2(1,i,j)], 1, file);
         fwrite(&U[c2(2,i,j)], sizeof U[c2(2,i,j)], 1, file);
         fwrite(&U[c2(3,i,j)], sizeof U[c2(3,i,j)], 1, file);
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

int Output3_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2, size_X3;
   char ext[20];
   char dat[20];
   char archivo[20];
   int num;
   
   size_X1 = Nx1-2*gc+1;
   size_X2 = Nx2-2*gc+1;
   size_X3 = Nx3-2*gc+1;
   
   num = *itprint;
   strcpy(ext,".bin");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   file = fopen(archivo,"wb");

   fwrite(&time, sizeof time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&size_X3, sizeof size_X3, 1, file);
   fwrite(&X1[gc], sizeof X1[gc], 1, file);
   fwrite(&X1[Nx1-gc], sizeof X1[gc], 1, file);   
   fwrite(&X2[gc], sizeof X2[gc], 1, file);
   fwrite(&X2[Nx2-gc], sizeof X2[gc], 1, file);  
   fwrite(&X3[gc], sizeof X3[gc], 1, file);
   fwrite(&X3[Nx3-gc], sizeof X3[gc], 1, file);  
         
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx2-gc; k++)
         {
            fwrite(&U[c3(0,i,j,k)], sizeof U[c3(0,i,j,k)], 1, file);
            fwrite(&U[c3(1,i,j,k)], sizeof U[c3(1,i,j,k)], 1, file);
            fwrite(&U[c3(2,i,j,k)], sizeof U[c3(2,i,j,k)], 1, file);
            fwrite(&U[c3(3,i,j,k)], sizeof U[c3(3,i,j,k)], 1, file);
            fwrite(&U[c3(4,i,j,k)], sizeof U[c3(4,i,j,k)], 1, file);
         }
      }
   }

   printf("itprint %d \n",*itprint);
   fclose(file);

   return 0;
}

