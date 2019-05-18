/**
 * @file restart.c
 *
 * @author Emilio Tejeda
 *
 * @brief Functions to restart from a given file
 */

//Do not erase any of these libraries//
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
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
#if dim == 1
   for(i = gc; i <= Nx1-gc; i++)
   {
      idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,\
      &U[c1(0,i)],&U[c1(1,i)],&U[c1(2,i)]);
   }

#elif dim == 2

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         idum = fscanf(file,"%lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
         &U[c2(0,i,j)],&U[c2(1,i,j)],&U[c2(2,i,j)],&U[c2(3,i,j)]);
      }
   }

#elif dim == 4

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         idum = fscanf(file,"%lf %lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
         &U[c2(0,i,j)],&U[c2(1,i,j)],&U[c2(2,i,j)],&U[c2(3,i,j)],&U[c2(4,i,j)]);
      }
   }

#elif dim == 3

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

#endif

   // CALLING BOUNDARIES TO GET GHOST CELLS RIGHT
    BOUNDARIES(U);

    fclose(file);
}                                        

void RESTART_BIN()
{
   FILE *file;
   int i, j, k, idum, ignore;
   double dum;
   char line[100];

   //Initialize dt
   dt = 0.1;

   file = fopen(restartfile,"rb");

   // Read time
   ignore = fread(&time, sizeof time, 1, file);

#if dim == 1

   // Read Nx1
   ignore = fread(&idum, sizeof idum, 1, file);

   if (idum + 2*gc - 1 != Nx1) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // check domain limits
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     

   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
        ignore = fread(&U[c1(0,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(1,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(2,i)], sizeof dum, 1, file);         
   }
    
#elif dim == 2 || dim == 4

   // Read Nx1
   ignore = fread(&idum, sizeof idum, 1, file);

   if (idum + 2*gc - 1 != Nx1) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // Read Nx2
   ignore = fread(&idum, sizeof idum, 1, file);
   
   if (idum + 2*gc - 1 != Nx2) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // check domain limits
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X2[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X2[Nx2-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }   
   
   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
        ignore = fread(&U[c1(0,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(1,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(2,i)], sizeof dum, 1, file);         
   }
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
        ignore = fread(&U[c2(0,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(1,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(2,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(3,i,j)], sizeof dum, 1, file);
      }
   }
    
#elif dim == 3

   // Read Nx1
   ignore = fread(&idum, sizeof idum, 1, file);

   if (idum + 2*gc - 1 != Nx1) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // Read Nx2
   ignore = fread(&idum, sizeof idum, 1, file);
   
   if (idum + 2*gc - 1 != Nx2) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // Read Nx3
   ignore = fread(&idum, sizeof idum, 1, file);
   
   if (idum + 2*gc - 1 != Nx3) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // check domain limits
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X2[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X2[Nx2-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }   
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X3[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != X3[Nx3-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }   
   
   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
        ignore = fread(&U[c1(0,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(1,i)], sizeof dum, 1, file);
        ignore = fread(&U[c1(2,i)], sizeof dum, 1, file);         
   }
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
        ignore = fread(&U[c2(0,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(1,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(2,i,j)], sizeof dum, 1, file);
        ignore = fread(&U[c2(3,i,j)], sizeof dum, 1, file);
      }
   }
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            ignore = fread(&U[c3(0,i,j,k)], sizeof dum, 1, file);
            ignore = fread(&U[c3(1,i,j,k)], sizeof dum, 1, file);
            ignore = fread(&U[c3(2,i,j,k)], sizeof dum, 1, file);
            ignore = fread(&U[c3(3,i,j,k)], sizeof dum, 1, file);               
            ignore = fread(&U[c3(4,i,j,k)], sizeof dum, 1, file);                 
         }
      }
   }

#endif
   
   // CALLING BOUNDARIES TO GET GHOST CELLS RIGHT
   BOUNDARIES(U);

   fclose(file);
}                                        
