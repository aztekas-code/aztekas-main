/**
 * @file messages.c
 * 
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Helpful functions for \a aztekas.
 */

#include"main.h"

void Check_Sim_Parameters()
{
#ifdef _OPENMP
   FILE *command;

   /*
    * Command for knowing the physical cores of the system
    */
   command = popen("grep '^core id' /proc/cpuinfo |sort -u|wc -l","r");
   fscanf(command,"%d",&MAX_NUM_THREADS);

   fclose(command);
#endif

   printf("\n");
   printf("aaaaa  zzzzz  ttttt  eeeee  k   k  aaaaa  sssss\n");
   printf("    a     zz    t    e   e  k  k       a  ss   \n");
   printf("aaaaa   zzz     t    eeeee  kkk    aaaaa  sssss\n");
   printf("a   a  zz       t    e      k  k   a   a     ss\n");
   printf("aaaaa  zzzzz    t    eeeee  k   k  aaaaa  sssss\n");
   printf("\n");
   printf("Running aztekas simulation...\n");
#ifdef _OPENMP
   printf("Parallel version operating with %d threads of %d available.\n",OMP_NUM,MAX_NUM_THREADS);
#else
   printf("Serial version.\n");
#endif
   printf("\n");

   // Print physics used
   if(PHYSICS == HD)  printf("Performing a HD simulation.\n");
   if(PHYSICS == RHD) printf("Performing a RHD simulation.\n");

   // Coordinates
#if PHYSICS == HD

   if(COORDINATES == CARTESIAN) printf("Cartesian grid (x,y,z).\n");
   if(COORDINATES == CYLINDRICAL) printf("Cylindrical grid (R,z,phi).\n");
   #if POLAR == FALSE
   if(COORDINATES == SPHERICAL) printf("Spherical grid (r,theta,phi).\n");
   #elif POLAR == TRUE
   if(COORDINATES == SPHERICAL) printf("Polar grid (R,phi).\n");
   #endif

#elif PHYSICS == RHD

   #if COORDINATES == CARTESIAN

   if(METRIC==USER) printf("Cartesian grid in a User defined space-time (x,y,z).\n");
   if(METRIC==MINK) printf("Cartesian grid in a Minkowski space-time (x,y,z).\n");
   if(METRIC==SCHW) printf("Cartesian grid in a Schwarzschild space-time (x,y,z).\n");
   if(METRIC==EF) printf("Cartesian grid in a Eddington-Finkelstein space-time (x,y,z).\n");
   if(METRIC==BL) printf("Cartesian grid in a Boyer-Lindquist space-time (x,y,z).\n");
   if(METRIC==KS) printf("Cartesian grid in a Kerr-Schild space-time (x,y,z).\n");

   #elif COORDINATES == CYLINDRICAL

   if(METRIC==USER) printf("Cylindrical grid in a User defined space-time (R,z,phi).\n");
   if(METRIC==MINK) printf("Cylindrical grid in a Minkowski space-time (R,z,phi).\n");
   if(METRIC==SCHW) printf("Cylindrical grid in a Schwarzschild space-time (R,z,phi).\n");
   if(METRIC==EF) printf("Cylindrical grid in a Eddington-Finkelstein space-time (R,z,phi).\n");
   if(METRIC==BL) printf("Cylindrical grid in a Boyer-Lindquist space-time (R,z,phi).\n");
   if(METRIC==KS) printf("Cylindrical grid in a Kerr-Schild space-time (R,z,phi).\n");

   #elif COORDINATES == SPHERICAL && POLAR == FALSE
   if(METRIC==USER) printf("Spherical grid in a User defined space-time (r,theta,phi).\n");
   if(METRIC==MINK) printf("Spherical grid in a Minkowski space-time (r,theta,phi).\n");
   if(METRIC==SCHW) printf("Spherical grid in a Schwarzschild space-time (r,theta,phi).\n");
   if(METRIC==EF) printf("Spherical grid in a Eddington-Finkelstein space-time (r,theta,phi).\n");
   if(METRIC==BL) printf("Spherical grid in a Boyer-Lindquist space-time (r,theta,phi).\n");
   if(METRIC==KS) printf("Spherical grid in a Kerr-Schild space-time (r,theta,phi).\n");
   #elif COORDINATES == SPHERICAL && POLAR == TRUE
   if(METRIC==USER) printf("Polar grid in a User defined space-time (R,phi).\n");
   if(METRIC==MINK) printf("Polar grid in a Minkowski space-time (R,phi).\n");
   if(METRIC==SCHW) printf("Polar grid in a Schwarzschild space-time (R,phi).\n");
   if(METRIC==EF) printf("Polar grid in a Eddington-Finkelstein space-time (R,phi).\n");
   if(METRIC==BL) printf("Polar grid in a Boyer-Lindquist space-time (R,phi).\n");
   if(METRIC==KS) printf("Polar grid in a Kerr-Schild space-time (R,phi).\n");
   #endif

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
   if(RECONST == GODUNOV)  printf("Using a zero-order piecewise reconstruction for the primitive variables.\n");
   if(RECONST == MINMOD)   printf("Using a first-order piecewise MINMOD reconstruction for the primitive variables.\n");
   if(RECONST == MC)       printf("Using a first-order piecewise MC reconstruction for the primitive variables.\n");
   if(RECONST == SUPERBEE) printf("Using a first-order piecewise SUPERBEE reconstruction for the primitive variables.\n");
   if(RECONST == WENO5)    printf("Using a fifth-order WENO5 reconstruction for the primitive variables.\n");

   // Flux solver
   if(FLUX == HLL)  printf("HLL Riemann solver.\n");
   if(FLUX == HLLC) printf("HLLC Riemann solver.\n");

   printf("\n");

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
#ifdef _OPENMP
   fprintf(file,"Parallel version operating with %d threads of %d available.\n",OMP_NUM,MAX_NUM_THREADS);
#else
   fprintf(file,"Serial version.\n");
#endif
   fprintf(file,"\n");

   // Print physics used
   if(PHYSICS == HD)  fprintf(file,"Performing a HD simulation.\n");
   if(PHYSICS == RHD) fprintf(file,"Performing a RHD simulation.\n");

   // Coordinates
#if PHYSICS == HD

   if(COORDINATES == CARTESIAN) fprintf(file,"Cartesian grid (x,y,z).\n");
   if(COORDINATES == CYLINDRICAL) fprintf(file,"Cylindrical grid (R,z,phi).\n");
   #if POLAR == FALSE
   if(COORDINATES == SPHERICAL) fprintf(file,"Spherical grid (r,theta,phi).\n");
   #elif POLAR == TRUE
   if(COORDINATES == SPHERICAL) fprintf(file,"Polar grid (R,phi).\n");
   #endif

#elif PHYSICS == RHD

   #if COORDINATES == CARTESIAN

   if(METRIC==USER) fprintf(file,"Cartesian grid in a User defined space-time (x,y,z).\n");
   if(METRIC==MINK) fprintf(file,"Cartesian grid in a Minkowski space-time (x,y,z).\n");
   if(METRIC==SCHW) fprintf(file,"Cartesian grid in a Schwarzschild space-time (x,y,z).\n");
   if(METRIC==EF) fprintf(file,"Cartesian grid in a Eddington-Finkelstein space-time (x,y,z).\n");
   if(METRIC==BL) fprintf(file,"Cartesian grid in a Boyer-Lindquist space-time (x,y,z).\n");
   if(METRIC==KS) fprintf(file,"Cartesian grid in a Kerr-Schild space-time (x,y,z).\n");

   #elif COORDINATES == CYLINDRICAL

   if(METRIC==USER) fprintf(file,"Cylindrical grid in a User defined space-time (R,z,phi).\n");
   if(METRIC==MINK) fprintf(file,"Cylindrical grid in a Minkowski space-time (R,z,phi).\n");
   if(METRIC==SCHW) fprintf(file,"Cylindrical grid in a Schwarzschild space-time (R,z,phi).\n");
   if(METRIC==EF) fprintf(file,"Cylindrical grid in a Eddington-Finkelstein space-time (R,z,phi).\n");
   if(METRIC==BL) fprintf(file,"Cylindrical grid in a Boyer-Lindquist space-time (R,z,phi).\n");
   if(METRIC==KS) fprintf(file,"Cylindrical grid in a Kerr-Schild space-time (R,z,phi).\n");

   #elif COORDINATES == SPHERICAL && POLAR == FALSE
   if(METRIC==USER) fprintf(file,"Spherical grid in a User defined space-time (r,theta,phi).\n");
   if(METRIC==MINK) fprintf(file,"Spherical grid in a Minkowski space-time (r,theta,phi).\n");
   if(METRIC==SCHW) fprintf(file,"Spherical grid in a Schwarzschild space-time (r,theta,phi).\n");
   if(METRIC==EF) fprintf(file,"Spherical grid in a Eddington-Finkelstein space-time (r,theta,phi).\n");
   if(METRIC==BL) fprintf(file,"Spherical grid in a Boyer-Lindquist space-time (r,theta,phi).\n");
   if(METRIC==KS) fprintf(file,"Spherical grid in a Kerr-Schild space-time (r,theta,phi).\n");
   #elif COORDINATES == SPHERICAL && POLAR == TRUE
   if(METRIC==USER) fprintf(file,"Polar grid in a User defined space-time (R,phi).\n");
   if(METRIC==MINK) fprintf(file,"Polar grid in a Minkowski space-time (R,phi).\n");
   if(METRIC==SCHW) fprintf(file,"Polar grid in a Schwarzschild space-time (R,phi).\n");
   if(METRIC==EF) fprintf(file,"Polar grid in a Eddington-Finkelstein space-time (R,phi).\n");
   if(METRIC==BL) fprintf(file,"Polar grid in a Boyer-Lindquist space-time (R,phi).\n");
   if(METRIC==KS) fprintf(file,"Polar grid in a Kerr-Schild space-time (R,phi).\n");
   #endif

#endif

   // Resolution
   if(DIM == 1) fprintf(file,"1D simulation with resolution %d grid cells\n",Nx1);
   if(DIM == 2) fprintf(file,"2D simulation with resolution %dX%d grid cells\n",Nx1,Nx2);
   if(DIM == 4) fprintf(file,"2.5D simulation with resolution %dX%d grid cells\n",Nx1,Nx2);
   if(DIM == 3) fprintf(file,"3D simulation with resolution %dX%dX%d grid cells\n",Nx1,Nx2,Nx3);

   // Equation of state
   if(EOS == IDEAL) fprintf(file,"Ideal equation of state with adiabatic index %f.\n",K);
   if(EOS == DUST)  fprintf(file,"Dust.\n");
   if(EOS == STIFF)  fprintf(file,"Stiff equation of state.\n");

   // Print MoL-RK order
   fprintf(file,"Time integration using a second order MoL-Runge Kutta.\n");

   // Print spatial numerical methods, algorithms and parameters.
   // Primitive variable reconstruction.
   if(RECONST == GODUNOV)  fprintf(file,"Using a zero-order piecewise reconstruction for the primitive variables.\n");
   if(RECONST == MINMOD)   fprintf(file,"Using a first-order piecewise MINMOD reconstruction for the primitive variables.\n");
   if(RECONST == MC)       fprintf(file,"Using a first-order piecewise MC reconstruction for the primitive variables.\n");
   if(RECONST == SUPERBEE) fprintf(file,"Using a first-order piecewise SUPERBEE reconstruction for the primitive variables.\n");
   if(RECONST == WENO5)    fprintf(file,"Using a fifth-order WENO5 reconstruction for the primitive variables.\n");

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

   #if MDOT == TRUE
   char create_mdot[] = "mkdir -p ";
   strcat(create_mdot,outputdirectory);
   strcat(create_mdot,"/Mdot/");
   sysret = system(create_mdot);
   #endif

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

void Check_Paramfile(char *param, int argc, char *argv[])
{
   if(argc != 2)
   {
      printf("%s\n","Wrong number of arguments") ;
      printf("%s\n","Execute as:") ;
      printf("%s\n","./aztekas paramfile") ;
      exit(EXIT_FAILURE);
   }

   strcpy(param, argv[1]);
}

void Computing_Time_Start()
{
#ifdef _OPENMP
   start = omp_get_wtime();
   omp_set_num_threads(OMP_NUM);
#else
   start = clock();
#endif
}

void Ending_Message()
{
   printf("\n");                                                                
   printf("AZTEKAS termination\n");                                             

#ifdef _OPENMP                                                                  
   if(omp_get_wtime() - start > 61)
   {
      int time_sec = (int)(omp_get_wtime()-start);                                 
      int hr  = time_sec/3600;                                                     
      int min = (time_sec%3600)/60;                                                
      int sec = (time_sec%60)%60;                                                  
      printf("Expend %d hr : %d min : %d sec in the parallelized version using %d threads of %d available.\n",hr,min,sec,OMP_NUM,MAX_NUM_THREADS);
      printf("\n");                                                                
   }
   else
   {
      printf("Expend %f sec in the parallelized version using %d threads of %d available.\n",omp_get_wtime()-start,OMP_NUM,MAX_NUM_THREADS);
      printf("\n");                                                                
   }
#else                                                                           
   if((double)(clock()-start)/CLOCKS_PER_SEC > 61.0)
   {
      int time_sec = (int)((double)(clock()-start)/CLOCKS_PER_SEC);                
      int hr  = time_sec/3600;                                                     
      int min = (time_sec%3600)/60;                                                
      int sec = (time_sec%60)%60;                                                  
      printf("Expend %d hr : %d min : %d sec in the serial version.\n",hr,min,sec);    
      printf("\n");                                                                
   }
   else
   {
      printf("Expend %f sec in the serial version.\n",(double)(clock()-start)/CLOCKS_PER_SEC);
      printf("\n");                                                                
   }
#endif                                                                          
}
