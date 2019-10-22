/**
 * @file auxfunc.c
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Helpful functions for \a aztekas.
 */

#include"main.h"

/**
 * Matrix times a vector
 */

int MxV(double *M, double *V, double *L)
{
   int n, m;
   double res=0.0;

   for(m = 0; m < eq; m++)
   {
      for(n = 0; n < eq; n++)
      {
         res += M[m*(eq) + n]*V[n];
      }

      L[m] = res;
      res = 0.0;
   }

   return 0;
}

void RoundGen(double *num)
{
   double r;
   double bla;
   double decnum;
   int expnum;
   
   if(*num != 0.0e+00)
   {
      expnum = (int)floor(log(fabs(*num))/log(10.0));
      decnum = *num/pow(10,(double)expnum);
      decnum = roundf(decnum*1.0e+15)/1.0e+15;
      *num = decnum*pow(10,(double)expnum);
   }
}

void Check_Sim_Parameters()
{
   MAX_NUM_THREADS = omp_get_max_threads();
   printf("\n");
   printf("aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n");
   printf("    a     zz    t    e   e  k  k       a  ss   \n");
   printf("aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n");
   printf("a   a  zz       t    e      k  k   a   a     ss\n");
   printf("aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n");
   printf("\n");
   printf("Running aztekas simulation...\n");
   printf("Operating with %d threads of %d available\n",OMP_NUM,MAX_NUM_THREADS);
   printf("\n");

   // Print physics used
   if(PHYSICS == HD)  printf("Performing a HD simulation.\n");
   if(PHYSICS == RHD) printf("Performing a RHD simulation.\n");

   // Coordinates

   #if PHYSICS == HD

   if(COORDINATES == CARTESIAN) printf("Cartesian grid.\n");
   if(COORDINATES == CYLINDRICAL) printf("Cylindrical grid.\n");
   if(COORDINATES == SPHERICAL) printf("Spherical grid.\n");

   #elif PHYSICS == RHD

   // Coordinates in RHD
   if(COORDINATES == CARTESIAN   && METRIC==USER) printf("Cartesian grid in a User defined space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==USER) printf("Cylindrical grid in a User defined space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==USER) printf("Spherical grid in a User defined space-time.\n");

   if(COORDINATES == CARTESIAN   && METRIC==MINK) printf("Cartesian grid in a Minkowski space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==MINK) printf("Cylindrical grid in a Minkowski space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==MINK) printf("Spherical grid in a Minkowski space-time.\n");

   if(COORDINATES == CARTESIAN   && METRIC==SCHW) printf("Cartesian grid in a Schwarzschild space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==SCHW) printf("Cylindrical grid in a Schwarzschild space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==SCHW) printf("Spherical grid in a Schwarzschild space-time.\n");

   if(COORDINATES == CARTESIAN   && METRIC==EF) printf("Cartesian grid in a Eddington-Finkelstein space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==EF) printf("Cylindrical grid in a Eddington-Finkelstein space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==EF) printf("Spherical grid in a Eddington-Finkelstein space-time.\n");

   if(COORDINATES == CARTESIAN   && METRIC==BL) printf("Cartesian grid in a Boyer-Lindquist space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==BL) printf("Cylindrical grid in a Boyer-Lindquist space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==BL) printf("Spherical grid in a Boyer-Lindquist space-time.\n");

   if(COORDINATES == CARTESIAN   && METRIC==KS) printf("Cartesian grid in a Kerr-Schild space-time.\n");
   if(COORDINATES == CYLINDRICAL && METRIC==KS) printf("Cylindrical grid in a Kerr-Schild space-time.\n");
   if(COORDINATES == SPHERICAL   && METRIC==KS) printf("Spherical grid in a Kerr-Schild space-time.\n");

   #endif

   // Resolution
   if(DIM == 1) printf("1D simulation with resolution %d grid cells\n",Nx1);
   if(DIM == 2) printf("2D simulation with resolution %dX%d grid cells\n",Nx1,Nx2);
   if(DIM == 4) printf("2.5D simulation with resolution %dX%d grid cells\n",Nx1,Nx2);
   if(DIM == 3) printf("3D simulation with resolution %dX%dX%d grid cells\n",Nx1,Nx2,Nx3);

   // Equation of state
   if(EOS == IDEAL) printf("Ideal equation of state with adiabatic index %f.\n",K);
   if(EOS == DUST)  printf("Dust.\n");
   if(EOS == STIFF)  printf("Stiff equation of state.\n");

   // Print MoL-RK order
   printf("Time integration using a second order MoL-Runge Kutta.\n");

   // Print spatial numerical methods, algorithms and parameters.
   // Primitive variable reconstruction.
   if(RECONST == GODUNOV)  printf("Zero-order piecewise reconstruction.\n");
   if(RECONST == MINMOD)   printf("MINMOD reconstruction.\n");
   if(RECONST == MC)       printf("MC reconstruction.\n");
   if(RECONST == SUPERBEE) printf("SUPERBEE reconstruction.\n");
   if(RECONST == WENO5)    printf("WENO5 reconstruction.\n");

   // Flux solver
   if(FLUX == HLL)  printf("HLL Riemann solver.\n");
   if(FLUX == HLLC) printf("HLLC Riemann solver.\n");

   FILE *file;

   char dum[100];
   strcat(dum,outputdirectory);
   strcat(dum,"/INFO");
   strcat(dum,"/info.sim");
   file = fopen(dum,"w");

   fprintf(file,"\n");
   fprintf(file,"aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n");
   fprintf(file,"    a     zz    t    e   e  k  k       a  ss   \n");
   fprintf(file,"aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n");
   fprintf(file,"a   a  zz       t    e      k  k   a   a     ss\n");
   fprintf(file,"aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n");
   fprintf(file,"\n");
   fprintf(file,"Running aztekas simulation...\n");
   fprintf(file,"\n");

   // Print physics used
   if(PHYSICS == HD)  fprintf(file,"Performing a HD simulation.\n");
   if(PHYSICS == RHD) fprintf(file,"Performing a RHD simulation.\n");

   // Coordinates
   if(COORDINATES == CARTESIAN) fprintf(file,"Cartesian grid.\n");
   if(COORDINATES == CYLINDRICAL) fprintf(file,"Cylindrical grid.\n");
   if(COORDINATES == SPHERICAL) fprintf(file,"Spherical grid.\n");

   // Resolution
   if(DIM == 1) fprintf(file,"1D simulation with resolution %d\n",Nx1);
   if(DIM == 2) fprintf(file,"2D simulation with resolution %dX%d\n",Nx1,Nx2);
   if(DIM == 4) fprintf(file,"2.5D simulation with resolution %dX%d\n",Nx1,Nx2);
   if(DIM == 3) fprintf(file,"3D simulation with resolution %dX%dX%d\n",Nx1,Nx2,Nx3);


   // Equation of state
   if(EOS == IDEAL) fprintf(file,"Ideal equation of state with adiabatic index %f.\n",K);
   if(EOS == DUST)  fprintf(file,"Dust.\n");
   if(EOS == STIFF) fprintf(file,"Stiff equation of state.\n");

   // Print MoL-RK order
   fprintf(file,"Time integration using a second order MoL-Runge Kutta.\n");

   // Print spatial numerical methods, algorithms and parameters.
   // Primitive variable reconstruction.
   if(RECONST == GODUNOV)  fprintf(file,"Zero-order piecewise reconstruction.\n");
   if(RECONST == MINMOD)   fprintf(file,"MINMOD reconstruction.\n");
   if(RECONST == MC)       fprintf(file,"MC reconstruction.\n");
   if(RECONST == SUPERBEE) fprintf(file,"SUPERBEE reconstruction.\n");
   if(RECONST == WENO5)    fprintf(file,"WENO5 reconstruction.\n");

   // Flux solver
   if(FLUX == HLL)  fprintf(file,"HLL Riemann solver.\n");
   if(FLUX == HLLC) fprintf(file,"HLLC Riemann solver.\n");

   fclose(file);
}

void Manage_Simulation_Info(int argc, char *argv[])
{
	// create output directory
   char create_dir[] = "mkdir -p ";	
	strcat(create_dir,outputdirectory);	
	int sysret = system(create_dir);

   char create_info[] = "mkdir -p ";
   strcat(create_info,outputdirectory);
   strcat(create_info,"/INFO/");
   sysret = system(create_info);

   // copy
   char dum[50] = "cp ";
   strcat(dum,paramfile_name);
   strcat(dum," ");
   strcat(dum,outputdirectory);
   strcat(dum,"INFO/");
   sysret = system(dum);
   
   strcpy(dum,"cp ");
   strcat(dum,"Makefile");
   strcat(dum," ");
   strcat(dum,outputdirectory);
   strcat(dum,"INFO/");
   sysret = system(dum);

   strcpy(dum,"cp ");
   strcat(dum,"*.c");
   strcat(dum," ");
   strcat(dum,outputdirectory);
   strcat(dum,"INFO/");
   sysret = system(dum);

   strcpy(dum,"cp ");
   strcat(dum,"*.h");
   strcat(dum," ");
   strcat(dum,outputdirectory);
   strcat(dum,"INFO/");
   sysret = system(dum);

   Check_Sim_Parameters();
   if(check_param == TRUE) getchar();
}
