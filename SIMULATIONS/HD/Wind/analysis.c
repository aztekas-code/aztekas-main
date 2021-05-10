/**
 * @file functions.c
 *
 * @brief 
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 22-04-2021 - 18:30:07
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 05-04-2021 - 16:20:36
 */
#include"main.h"

void Analysis(double *B)
{
   double P[3];
   FILE *file;
   char archivo[50];
   char dat[50];
   eos_ eos;
   gauge_ local_grid;

   snprintf(dat,8,"%d",itprint);
   strcpy(archivo,outputdirectory);
   strcat(archivo,"analysis_");
   strcat(archivo,dat);
   strcat(archivo,".dat");

   if(grid.time >= tprint)
   {
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
                                                                                   
      fprintf(file,"#X1, X2, TEMPERATURE, SPEED OF SOUND, ENERGY, ENTROPY\n");  

      for(int i = gc; i <= Nx1-gc; i++)
      {
         for(int j = gc; j <= Nx2-gc; j++)
         {
            P[0] = B(RHO,i,j);
            P[1] = B(PRE,i,j);
            P[2] = 0.0;

            EoS(&eos,P,&local_grid);

            fprintf(file,"%e %e %e %e %e %e\n", grid.X1[i],grid.X2[j],B(RHO,i,j),B(PRE,i,j),eos.temp,eos.cs);
         }
      }

      fclose(file);
   }
}
