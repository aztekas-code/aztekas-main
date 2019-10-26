/**
 * @file print_values.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Output functions: ASCII and Binary.
 */

#include"main.h"

#if DIM == 1

void Output1(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[50];
   char dat[50];
   char archivo[50];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"w");   

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",grid.time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
#if COORDINATES == CARTESIAN
   fprintf(file,"CARTESIAN\n");
#elif COORDINATES == CYLINDRICAL
   fprintf(file,"CYLINDRICAL\n");
#elif COORDINATES == SPHERICAL
   fprintf(file,"SPHERICAL \n");
#elif POLAR == TRUE
   fprintf(file,"POLAR \n");
#endif
   fprintf(file,"###################################\n");

   for(i = gc; i <= Nx1-gc; i++)
   {
      fprintf(file,"%e %e %e %e\n",grid.X1[i],\
      U(0,i),\
      U(1,i),\
      U(2,i));
   }

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

void Output1_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1;
   char ext[50];
   char dat[50];
   char archivo[50];
   int num;
   
   size_X1 = Nx1-2*gc+1;
   
   num = *itprint;
   strcpy(ext,".bin");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"wb");

   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
         
   for(i = gc; i <= Nx1-gc; i++)
	{
      fwrite(&U(0,i), sizeof U(0,i), 1, file);
      fwrite(&U(1,i), sizeof U(1,i), 1, file);
      fwrite(&U(2,i), sizeof U(2,i), 1, file);
   }

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

#elif DIM == 2 || DIM == 4

void Output2(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[50];
   char dat[50];
   char archivo[50];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"w");

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",grid.time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
#if COORDINATES == CARTESIAN
   fprintf(file,"CARTESIAN\n");
#elif COORDINATES == CYLINDRICAL
   fprintf(file,"CYLINDRICAL\n");
#elif COORDINATES == SPHERICAL
   fprintf(file,"SPHERICAL \n");
#elif POLAR == TRUE
   fprintf(file,"POLAR \n");
#endif
   fprintf(file,"###################################\n");

#if DIM == 2
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fprintf(file,"%e %e %e %e %e %e\n",grid.X1[i],grid.X2[j],U(0,i,j),\
         U(1,i,j),\
         U(2,i,j),\
         U(3,i,j));
      }
   }
#elif DIM == 4
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fprintf(file,"%e %e %e %e %e %e %e\n",grid.X1[i],grid.X2[j],U(0,i,j),\
         U(1,i,j),\
         U(2,i,j),\
         U(3,i,j),\
         U(4,i,j));
      }
   }
#endif 

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

void Output2_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2;
   char ext[50];
   char dat[50];
   char archivo[50];
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
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"wb");

   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
   fwrite(&grid.X2[gc], sizeof grid.X2[gc], 1, file);
   fwrite(&grid.X2[Nx2-gc], sizeof grid.X2[gc], 1, file);  
         
#if DIM == 2
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fwrite(&U(0,i,j), sizeof U(0,i,j), 1, file);
         fwrite(&U(1,i,j), sizeof U(1,i,j), 1, file);
         fwrite(&U(2,i,j), sizeof U(2,i,j), 1, file);
         fwrite(&U(3,i,j), sizeof U(3,i,j), 1, file);
      }
   }
#elif DIM == 4
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fwrite(&U(0,i,j), sizeof U(0,i,j), 1, file);
         fwrite(&U(1,i,j), sizeof U(1,i,j), 1, file);
         fwrite(&U(2,i,j), sizeof U(2,i,j), 1, file);
         fwrite(&U(3,i,j), sizeof U(3,i,j), 1, file);
         fwrite(&U(4,i,j), sizeof U(4,i,j), 1, file);
      }
   }
#endif

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

#elif DIM == 3

void Output3(int *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[50];
   char dat[50];
   char archivo[50];
   int num;

   num = *itprint;
   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);
   strcpy(archivo,outputdirectory);
   strcat(archivo,outputfile);
   strcat(archivo,dat);
   strcat(archivo,ext);
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"w");

   fprintf(file,"###############PARAM###############\n");
   fprintf(file,"%e \n",grid.time);
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
   fprintf(file,"%d \n",Nx3-2*gc+1);
#if COORDINATES == CARTESIAN
   fprintf(file,"CARTESIAN\n");
#elif COORDINATES == CYLINDRICAL
   fprintf(file,"CYLINDRICAL\n");
#elif COORDINATES == SPHERICAL
   fprintf(file,"SPHERICAL \n");
#elif POLAR == TRUE
   fprintf(file,"POLAR \n");
#endif
   fprintf(file,"###################################\n");

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = Nx2/2; j <= Nx2/2; j++)
      {
         for(k = Nx3/2; k <= Nx3/2; k++)
         {
            fprintf(file,"%f %f %f %f %f %f %f %f\n",grid.X1[i],grid.X2[j],grid.X3[k],\
            U(0,i,j,k),\
            U(1,i,j,k),\
            U(2,i,j,k),\
            U(3,i,j,k),\
            U(4,i,j,k));
         }
      }
   }

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

void Output3_bin(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2, size_X3;
   char ext[50];
   char dat[50];
   char archivo[50];
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
   if(CHECK_NAN == TRUE)
   {
      strcpy(archivo,outputdirectory);
      strcat(archivo,"last.dat");   
   }
   file = fopen(archivo,"wb");

   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&size_X3, sizeof size_X3, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
   fwrite(&grid.X2[gc], sizeof grid.X2[gc], 1, file);
   fwrite(&grid.X2[Nx2-gc], sizeof grid.X2[gc], 1, file);  
   fwrite(&grid.X3[gc], sizeof grid.X3[gc], 1, file);
   fwrite(&grid.X3[Nx3-gc], sizeof grid.X3[gc], 1, file);  
         
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx2-gc; k++)
         {
            fwrite(&U(0,i,j,k), sizeof U(0,i,j,k), 1, file);
            fwrite(&U(1,i,j,k), sizeof U(1,i,j,k), 1, file);
            fwrite(&U(2,i,j,k), sizeof U(2,i,j,k), 1, file);
            fwrite(&U(3,i,j,k), sizeof U(3,i,j,k), 1, file);
            fwrite(&U(4,i,j,k), sizeof U(4,i,j,k), 1, file);
         }
      }
   }

   printf("itprint : %d, output file : %s\n",*itprint,archivo);
   fclose(file);
}

#endif 
