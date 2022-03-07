/*
 * aztekas initial module
 * Date of creation/modification: 27-09-19 11:26:00
 * author: Alejandro Aguayo-Ortiz
 */

//Do not erase any of these libraries//
#include"main.h"

void Initial()
{
   FILE *file;
   char file_name[50];
   double u[3];
   eos_ eos;
   gauge_ local_grid;

   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   ///////////////////////////
   //---------Wind----------//
   ///////////////////////////

#if EOS == HELMHOLTZ
   u[0] = density_inf;
   u[1] = temperature_inf;
   u[2] = 0.0;

   EoS_DT(&eos,u);

//   printf("Dens %.16lf\n",eos.rho);
//   printf("Pre %.16lf\n",eos.p);
//   printf("Temp %.16lf\n",eos.temp);
//   printf("Cs %e\n",eos.cs);
//   getchar();

   pressure_inf = eos.p;
   Mach_inf     = velocity_inf/eos.cs;

//   printf("Speed of sound c_s = %e\n",eos.cs*vel_units);
//   printf("for T = %e\n",temperature_inf);
//   getchar();
#elif EOS == IDEAL

   pressure_inf = (density_inf * dens_units) * K_B_cgs * (temperature_inf * temp_units) / (dens_units * vel_units * vel_units * mH_cgs);
   eos.cs     = sqrt(K * pressure_inf / density_inf);

//   Mach_inf = velocity_inf/eos.cs;
   velocity_inf = Mach_inf * eos.cs;

#endif

   velocity_inf = sqrt(2.0);
   eos.cs = velocity_inf / Mach_inf;
   pressure_inf = density_inf * eos.cs * eos.cs / K;

   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) =  density_inf;
         U(PRE,i,j) =  pressure_inf;
         U(VX1,i,j) =  velocity_inf*cos(grid.X2[j]);
         U(VX2,i,j) = -velocity_inf*sin(grid.X2[j]);
      }
   }

   strcpy(file_name,outputdirectory);
   strcat(file_name,"info.sim");

   printf("Printing CGS info in %s\n",file_name);
   printf("\n");

   file = fopen(file_name,"w");

   fprintf(file,"Density at boundary = %e g cm⁻³\n",density_inf*dens_units);
   fprintf(file,"Pressure at boundary = %e erg cm⁻³\n",pressure_inf*dens_units*vel_units*vel_units);
   fprintf(file,"Speed of sound at boundary = %e cm s⁻¹\n",eos.cs*vel_units);
   fprintf(file,"Wind velocity = %e cm s⁻¹\n",velocity_inf*vel_units);
   fprintf(file,"Mach number of wind = %e\n",Mach_inf);
   fprintf(file,"\n");
   fprintf(file,"Density in units of: %e g cm⁻³\n",dens_units);
   fprintf(file,"Pressure in units of: %e erg cm⁻³\n",dens_units*vel_units*vel_units);
   fprintf(file,"Velocity in units of: %e g cm⁻³\n",vel_units);
   fprintf(file,"Temperature in units of: %e K\n",temp_units);
   fprintf(file,"Distances in units of 2GM/%e cm\n",vel_units*vel_units);
   fprintf(file,"Time in units of GM/%e sec\n",vel_units*vel_units*vel_units);

   fclose(file);

   double density, P_goal, t_guess;
   double delta_P, delta_T, t_lower, t_upper;
   double err_P, err_T;
   int count = 1;
   int max_iter = 100;

   //density = 1.000000e+00;
   //P_goal  = 8.287063e-02;
//   density = 1.0;
//   P_goal  = 1.0;
//   t_guess = 1.000000e+00;
//
//   t_lower = 1.0e-03;
//   t_upper = 1.0e+05;
//
//   printf("Density %e\n",density);
//   printf("P_goal  %e\n",P_goal);
//   printf("T_guess %e\n",t_guess);
//
//   while (count < max_iter)
//   {
//      u[0] = density;
//      u[1] = t_guess;
//
//      EoS_DT(&eos,u);
//
//      delta_P = P_goal - eos.p;
//
//      if(delta_P > 0.0)
//      {
//         t_lower = eos.temp;
//      }
//      else
//      {
//         t_upper = eos.temp;
//      }
//
//      delta_T  = delta_P / eos.dpt;
//      eos.temp = eos.temp + delta_T;
//
//      if (eos.temp > t_upper || eos.temp < t_lower)
//         eos.temp = sqrt(t_upper*t_lower);
//
//      err_P = delta_P / P_goal;
//      err_T = delta_T / eos.temp;
//
//      count++;
//      printf("%e %e %e %e\n",t_guess,delta_T,delta_P,eos.dpt);
//      printf("%e %e\n",err_T,err_P);
//      printf("%e %e %e\n",eos.temp,eos.rho,eos.p);
//      printf("%d\n",count);
//
//      if((fabs(err_P) <= 1.0e-04 && fabs(err_T) <= 1.0e-04) || fabs(1.0 - eos.temp/t_guess) < 1.0e-04)
//         count = max_iter;
//
//      t_guess = eos.temp;
//
//      getchar();
//   }
}
