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

   pressure_inf = eos.p;
   Mach_inf     = velocity_inf/eos.cs;

//   printf("Speed of sound c_s = %e\n",eos.cs*vel_units);
//   printf("for T = %e\n",temperature_inf);
//   getchar();
#elif EOS == IDEAL

   pressure_inf = (density_inf * dens_units) * K_B_cgs * (temperature_inf * temp_units) / (dens_units * vel_units * vel_units * mH_cgs);
   eos.cs     = sqrt(K * pressure_inf / density_inf);

   Mach_inf = velocity_inf/eos.cs;

#endif

   for(int i = 0; i <= Nx1; i++)
   {
      for(int j = 0; j <= Nx2; j++)
      {
         U(RHO,i,j) =  density_inf;
         U(PRE,i,j) =  pressure_inf;//U(RHO,i,j)/K;//pressure_inf;//pow(U(RHO,i,j),K)/K;
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
   fprintf(file,"Distances in units of GM/%e cm\n",vel_units*vel_units);
   fprintf(file,"Time in units of GM/%e sec\n",vel_units*vel_units*vel_units);

   fclose(file);
}
