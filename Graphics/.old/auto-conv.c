#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int main(int argc, char* argv[])
{
   FILE *fdata;
   char line[100];
   char file[100];
   int i, j;
   int Nx1, Nx2;
   int idum;
   double time, dum;

   if(argc != 3)
   {  
      printf("%s\n", "Wrong number of arguments");
      printf("%s\n", "./acc file file file");
      exit(EXIT_FAILURE);
   }

   strcpy(file,argv[1]);

   fdata = fopen(file,"r");
   idum = fscanf(fdata,"%s\n",line);
   idum = fscanf(fdata,"%lf\n",&time); 
   idum = fscanf(fdata,"%d\n",&Nx1);
   idum = fscanf(fdata,"%d\n",&Nx2);
   idum = fscanf(fdata,"%s\n",line);

   double *Radius1, *Theta1, *Density1;
   Radius1  = (double *)malloc(Nx1*Nx2*sizeof(double));
   Theta1   = (double *)malloc(Nx1*Nx2*sizeof(double));
   Density1 = (double *)malloc(Nx1*Nx2*sizeof(double));
   
   Nx1 = Nx1-1;
   Nx2 = Nx2-1;

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         idum = fscanf(fdata,"%lf %lf %lf %lf %lf %lf\n",&Radius1[i*(Nx2+1) + j],\
         &Theta1[i*(Nx2+1) + j],\
         &Density1[i*(Nx2+1) + j],\
         &dum,\
         &dum,\
         &dum);
      }
   }

   fclose(fdata);

   double *Density11, *Radius11, *Theta11;
   Radius11 = (double *)malloc(100000*sizeof(double));
   Theta11 = (double *)malloc(100000*sizeof(double));
   Density11 = (double *)malloc(100000*sizeof(double));

   for(i = 0; i <= Nx1; i++)
   {
      Radius11[i] = Radius1[(i)*(Nx2+1) + (Nx2)];
      Theta11[i] = Theta1[(i)*(Nx2+1) + (Nx2)];
      Density11[i] = Density1[(i)*(Nx2+1) + (Nx2)];
   }

   strcpy(file,argv[2]);

   fdata = fopen(file,"r");
   idum = fscanf(fdata,"%s\n",line);
   idum = fscanf(fdata,"%lf\n",&time); 
   idum = fscanf(fdata,"%d\n",&Nx1);
   idum = fscanf(fdata,"%d\n",&Nx2);
   idum = fscanf(fdata,"%s\n",line);

   double *Radius2, *Theta2, *Density2;
   Radius2  = (double *)malloc(Nx1*Nx2*sizeof(double));
   Theta2   = (double *)malloc(Nx1*Nx2*sizeof(double));
   Density2 = (double *)malloc(Nx1*Nx2*sizeof(double));

   Nx1 = (Nx1-1);
   Nx2 = (Nx2-1);

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         idum = fscanf(fdata,"%lf %lf %lf %lf %lf %lf\n",&Radius2[(i)*(Nx2+1) + (j)],\
         &Theta2[(i)*(Nx2+1) + (j)],\
         &Density2[(i)*(Nx2+1) + (j)],\
         &dum,\
         &dum,\
         &dum);
      }
   }

   fclose(fdata);
   
   double *Density22, *Radius22, *Theta22;
   Radius22 = (double *)malloc(100000*sizeof(double));
   Theta22 = (double *)malloc(100000*sizeof(double));
   Density22 = (double *)malloc(100000*sizeof(double));

   for(i = 0; i <= Nx1; i=i+2)
   {
      Density22[i/2] = Density2[(i)*(Nx2+1) + (Nx2)];
      Radius22[i/2] = Radius2[(i)*(Nx2+1) + (Nx2)];
      Theta22[i/2] = Theta2[(i)*(Nx2+1) + (Nx2)];
   }

   double sum = 0;

   for(i = 0; i <= Nx1/2; i++)
   {
      sum = sum + fabs(Density22[i] - Density11[i]);
      //printf("%e %e %e %e %e %e\n",Radius11[i],Radius22[i],Theta11[i],Theta22[i],Density11[i],Density22[i]);
   }

   printf("%e %e \n",time,2*sum*(38.46 - 1.9)/Nx1/1.0e-10);
   
   return 0;
}
