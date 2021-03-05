/**
 * @file restart.c
 *
 * @author Emilio Tejeda
 *
 * @brief Functions to restart from a given file
 */

#include"main.h"

void Restart()
{
   FILE *file;
   int i, j, k, idum, size;
   double dum;
   char line[100];

   //Initialize dt
   dt = 0.1;

   file = fopen(restartfile,"r");

   // Read rest of file an initialize variables      
#if DIM == 1
   // Skip first line
   idum = fscanf(file,"%s\n",line) ;

   // Read grid.time
   idum = fscanf(file,"%lf\n",&grid.time) ;

   // Read Nx1
   idum = fscanf(file,"%d\n",&size) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   for(i = gc; i <= Nx1-gc; i++)
   {
      idum = fscanf(file,"%lf %lf %lf %lf\n",&dum,\
      &U(0,i),&U(1,i),&U(2,i));
   }

#elif DIM == 2
   // Skip first line
   idum = fscanf(file,"%s\n",line) ;

   // Read grid.time
   idum = fscanf(file,"%lf\n",&grid.time) ;

   // Read Nx1
   idum = fscanf(file,"%d\n",&size) ;

   // Read Nx2
   idum = fscanf(file,"%d\n",&size) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         idum = fscanf(file,"%lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
         &U(RHO,i,j),&U(PRE,i,j),&U(VX1,i,j),&U(VX2,i,j));
      }
   }

#elif DIM == 4
   // Skip first line
   idum = fscanf(file,"%s\n",line) ;

   // Read grid.time
   idum = fscanf(file,"%lf\n",&grid.time) ;

   // Read Nx1
   idum = fscanf(file,"%d\n",&size) ;

   // Read Nx2
   idum = fscanf(file,"%d\n",&size) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   // Skip line   
   idum = fscanf(file,"%s\n",line) ;

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         idum = fscanf(file,"%lf %lf %lf %lf %lf %lf %lf\n",&dum,&dum,\
         &U(0,i,j),&U(1,i,j),&U(2,i,j),&U(3,i,j),&U(4,i,j));
      }
   }

#elif DIM == 3

   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            idum = fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf\n",\
            &dum,&dum,&dum,\
            &U(0,i,j,k),&U(1,i,j,k),\
            &U(2,i,j,k),&U(3,i,j,k),&U(4,i,j,k));
         }
      }
   }

#endif

   fclose(file);

   if(restart_filecount == 0)
   {
      grid.time = 0.0;
   }
}                                        

void Restart_Bin()
{
   FILE *file;
   int i, j, k, idum, ignore;
   double dum;
   char line[100];

   //Initialize dt
   dt = 0.1;

   file = fopen(restartfile,"rb");

   // Read grid.time
   ignore = fread(&grid.time, sizeof grid.time, 1, file);

#if DIM == 1

   // Read Nx1
   ignore = fread(&idum, sizeof idum, 1, file);

   if (idum + 2*gc - 1 != Nx1) 
   {
      printf("Size of input file incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }

   // check domain limits
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     

   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
        ignore = fread(&U(0,i), sizeof dum, 1, file);
        ignore = fread(&U(1,i), sizeof dum, 1, file);
        ignore = fread(&U(2,i), sizeof dum, 1, file);         
   }
    
#elif DIM == 2 

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
   if (dum != grid.X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[Nx2-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }   
   
   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
        ignore = fread(&U(0,i,j), sizeof dum, 1, file);
        ignore = fread(&U(1,i,j), sizeof dum, 1, file);
        ignore = fread(&U(2,i,j), sizeof dum, 1, file);
        ignore = fread(&U(3,i,j), sizeof dum, 1, file);
      }
   }
    
#elif DIM == 4 

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
   if (dum != grid.X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[Nx2-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }   
   
   // Read rest of file and initialize variables      
   for(i = gc; i <= Nx1-gc; i++)
   {
      for(j = gc; j <= Nx2-gc; j++)
      {
        ignore = fread(&U(0,i,j), sizeof dum, 1, file);
        ignore = fread(&U(1,i,j), sizeof dum, 1, file);
        ignore = fread(&U(2,i,j), sizeof dum, 1, file);
        ignore = fread(&U(3,i,j), sizeof dum, 1, file);
        ignore = fread(&U(4,i,j), sizeof dum, 1, file);
      }
   }
    
#elif DIM == 3

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
   if (dum != grid.X1[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X1[Nx1-gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[gc]) 
   {
      printf("Domain size incompatible with simulation parameters!\n");
      exit(EXIT_FAILURE);
   }     
   ignore = fread(&dum, sizeof dum, 1, file);
   if (dum != grid.X2[Nx2-gc]) 
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
      for(j = gc; j <= Nx2-gc; j++)
      {
         for(k = gc; k <= Nx3-gc; k++)
         {
            ignore = fread(&U(0,i,j,k), sizeof dum, 1, file);
            ignore = fread(&U(1,i,j,k), sizeof dum, 1, file);
            ignore = fread(&U(2,i,j,k), sizeof dum, 1, file);
            ignore = fread(&U(3,i,j,k), sizeof dum, 1, file);               
            ignore = fread(&U(4,i,j,k), sizeof dum, 1, file);                 
         }
      }
   }

#endif

   fclose(file);
}                                        
