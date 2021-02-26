/**
 * @file output.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Output functions: ASCII and Binary.
 */

#include"main.h"

void Output_ascii_int(int *itprint)
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

#if DIM == 1
   fprintf(file,"%d \n",Nx1-2*gc+1);
#elif DIM == 2 || DIM == 4
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
#elif DIM == 3
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
   fprintf(file,"%d \n",Nx3-2*gc+1);
#endif

#if COORDINATES == CARTESIAN
   fprintf(file,"CARTESIAN\n");
#elif COORDINATES == CYLINDRICAL
   fprintf(file,"CYLINDRICAL\n");
#elif COORDINATES == SPHERICAL && POLAR == FALSE
   fprintf(file,"SPHERICAL\n");
#elif POLAR == TRUE
   fprintf(file,"POLAR\n");
#endif

   fprintf(file,"###################################\n");

#if DIM == 1
   for(i = gc; i <= Nx1-gc; i++)
   {
      fprintf(file,"%e ",grid.X1[i]);
      for(n = 0; n < eq; n++)
      {
         if(n < eq - 1)
         {
            fprintf(file,"%e ",U(n,i));
         }
         if(n == eq - 1)
         {
            fprintf(file,"%e\n",U(n,i));
         }
      }
   }
#elif DIM == 2 || DIM == 4
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fprintf(file,"%e %e ",grid.X1[i],grid.X2[j]);
         for(n = 0; n < eq; n++)
         {
            if(n < eq - 1)
            {
               fprintf(file,"%e ",U(n,i,j));
            }
            if(n == eq - 1)
            {
               fprintf(file,"%e\n",U(n,i,j));
            }
         }
      }
   }
#elif DIM == 3
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            fprintf(file,"%e %e %e ",grid.X1[i],grid.X2[j],grid.X3[k]);
            for(n = 0; n < eq; n++)
            {
               if(n < eq - 1)
               {
                  fprintf(file,"%e ",U(n,i,j,k));
               }
               if(n == eq - 1)
               {
                  fprintf(file,"%e\n",U(n,i,j,k));
               }
            }
         }
      }
   }
#endif

   printf("id : %d, output file : %s\n",itprint,archivo);
   #if FILE_NOTIFICATION == TRUE
      char *a = "notification ";
      char script[100];
      strcpy(script,a);
      strcat(script,archivo);
      system(script);
   #endif

   fclose(file);
}

void Output_bin_int(int *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2, size_X3;
   char ext[50];
   char dat[50];
   char archivo[50];
   int num;
   
#if DIM == 1
   size_X1 = Nx1-2*gc+1;
#elif DIM == 2 || DIM == 4
   size_X1 = Nx1-2*gc+1;
   size_X2 = Nx2-2*gc+1;
#endif
   
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

#if DIM == 1
   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
         
   for(i = gc; i <= Nx1-gc; i++)
	{
      for(n = 0; n < eq; n++)
         fwrite(&U(n,i), sizeof U(n,i), 1, file);
   }
#elif DIM == 2 || DIM == 4
   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
   fwrite(&grid.X2[gc], sizeof grid.X2[gc], 1, file);
   fwrite(&grid.X2[Nx2-gc], sizeof grid.X2[gc], 1, file);  
         
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(n = 0; n < eq; n++)
            fwrite(&U(n,i,j), sizeof U(n,i,j), 1, file);
      }
   }
#endif

   printf("id : %d, output file : %s\n",itprint,archivo);
   #if FILE_NOTIFICATION == TRUE
      char *a = "notification ";
      char script[100];
      strcpy(script,a);
      strcat(script,archivo);
      system(script);
   #endif

   fclose(file);
}

void Output_ascii_char(char *itprint)
{
   FILE *file;
   int n, i, j, k;
   char ext[50];
   char dat[50];
   char archivo[50];

   strcpy(ext,".dat");
   strcpy(dat,itprint);
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

#if DIM == 1
   fprintf(file,"%d \n",Nx1-2*gc+1);
#elif DIM == 2 || DIM == 4
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
#elif DIM == 3
   fprintf(file,"%d \n",Nx1-2*gc+1);
   fprintf(file,"%d \n",Nx2-2*gc+1);
   fprintf(file,"%d \n",Nx3-2*gc+1);
#endif

#if COORDINATES == CARTESIAN
   fprintf(file,"CARTESIAN\n");
#elif COORDINATES == CYLINDRICAL
   fprintf(file,"CYLINDRICAL\n");
#elif COORDINATES == SPHERICAL && POLAR == FALSE
   fprintf(file,"SPHERICAL\n");
#elif POLAR == TRUE
   fprintf(file,"POLAR\n");
#endif

   fprintf(file,"###################################\n");

#if DIM == 1
   for(i = gc; i <= Nx1-gc; i++)
   {
      fprintf(file,"%e ",grid.X1[i]);
      for(n = 0; n < eq; n++)
      {
         if(n < eq - 1)
         {
            fprintf(file,"%e ",U(n,i));
         }
         if(n == eq - 1)
         {
            fprintf(file,"%e\n",U(n,i));
         }
      }
   }
#elif DIM == 2 || DIM == 4
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         fprintf(file,"%e %e ",grid.X1[i],grid.X2[j]);
         for(n = 0; n < eq; n++)
         {
            if(n < eq - 1)
            {
               fprintf(file,"%e ",U(n,i,j));
            }
            if(n == eq - 1)
            {
               fprintf(file,"%e\n",U(n,i,j));
            }
         }
      }
   }
#elif DIM == 3
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            fprintf(file,"%e %e %e ",grid.X1[i],grid.X2[j],grid.X3[k]);
            for(n = 0; n < eq; n++)
            {
               if(n < eq - 1)
               {
                  fprintf(file,"%e ",U(n,i,j,k));
               }
               if(n == eq - 1)
               {
                  fprintf(file,"%e\n",U(n,i,j,k));
               }
            }
         }
      }
   }
#endif

   printf("id : %s, output file : %s\n",itprint,archivo);
   #if FILE_NOTIFICATION == TRUE
      char *a = "notification ";
      char script[100];
      strcpy(script,a);
      strcat(script,archivo);
      system(script);
   #endif

   fclose(file);
}

void Output_bin_char(char *itprint)
{
   FILE *file;
   int n, i, j, k, size_X1, size_X2, size_X3;
   char ext[50];
   char dat[50];
   char archivo[50];
   
#if DIM == 1
   size_X1 = Nx1-2*gc+1;
#elif DIM == 2 || DIM == 4
   size_X1 = Nx1-2*gc+1;
   size_X2 = Nx2-2*gc+1;
#endif
   
   strcpy(ext,".bin");
   strcpy(dat,itprint);
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

#if DIM == 1
   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
         
   for(i = gc; i <= Nx1-gc; i++)
	{
      for(n = 0; n < eq; n++)
         fwrite(&U(n,i), sizeof U(n,i), 1, file);
   }
#elif DIM == 2 || DIM == 4
   fwrite(&grid.time, sizeof grid.time, 1, file);
   fwrite(&size_X1, sizeof size_X1, 1, file);
   fwrite(&size_X2, sizeof size_X2, 1, file);
   fwrite(&grid.X1[gc], sizeof grid.X1[gc], 1, file);
   fwrite(&grid.X1[Nx1-gc], sizeof grid.X1[gc], 1, file);   
   fwrite(&grid.X2[gc], sizeof grid.X2[gc], 1, file);
   fwrite(&grid.X2[Nx2-gc], sizeof grid.X2[gc], 1, file);  
         
   for(i = gc; i <= Nx1-gc; i++)
	 {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(n = 0; n < eq; n++)
            fwrite(&U(n,i,j), sizeof U(n,i,j), 1, file);
      }
   }
#endif

   printf("id : %s, output file : %s\n",itprint,archivo);
   #if FILE_NOTIFICATION == TRUE
      char *a = "notification ";
      char script[100];
      strcpy(script,a);
      strcat(script,archivo);
      system(script);
   #endif

   fclose(file);
}
