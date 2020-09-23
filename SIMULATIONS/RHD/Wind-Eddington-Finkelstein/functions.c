/**
 * @file functions.c
 * 
 * @brief Useful functions
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 22-09-2020 - 21:19:07
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 28-04-2020
 */

#include"main.h"

void Mass_Accretion_Rate(double *B)
{
   FILE *file_time, *file_domain;
   char name_time[100], name_domain[100];
   char ext[50], dat[50];
   int num;

   gauge_ local_grid;
   #if COORDINATES == SPHERICAL
   double time, r, th, thm, phi;
   double y_rr, y_thth, y_phiphi;
   double y_rth, y_rphi, y_thphi;
   double beta_r, beta_th, beta_phi;
   double lapse, dety;
   double rho, pre, vr, vth, vphi;
   double Vr, Vth, Vphi;
   double VV, W, Ur;
   double Mdot, Min, Mej, Vol;
   double Mean[Nx1+1];
   double Mdot_Acc = 0.0;
   double Mdot_Mean_Time = 0.0;
   double Delta;
   #endif

   local_grid.x[0] = grid.time;
   time            = local_grid.x[0];
   Mdot = 0.0;

   num = itprint;

   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);

   strcpy(name_time,outputdirectory);
   strcat(name_time,"Mdot/");
   strcat(name_time,"Mdot_Mean_Time");
   strcat(name_time,ext);

   strcpy(name_domain,outputdirectory);
   strcat(name_domain,"Mdot/");
   strcat(name_domain,"Mdot_Mean_Domain_");
   strcat(name_domain,dat);
   strcat(name_domain,ext);

   Mdot = 0.0;
   Min  = 0.0;
   Mej  = 0.0;

   for(int i = gc; i <= Nx1-gc; i++)
   {
      for(int j = gc+1; j <= Nx2-gc-1; j++)
      {
      #if DIM == 1                                                                    
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = 0.0;                                                       
         local_grid.x[3] = 0.0;                                                       
         #if COORDINATES == SPHERICAL                                                 
         local_grid.x[2] = M_PI_2;                                                    
         #endif                                                                       
                                                                                      
      #elif DIM == 2 || DIM == 4                                                      
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = grid.X2[j];                                                
         local_grid.x[3] = 0.0;                                                       
         #if POLAR == TRUE                                                            
         local_grid.x[2] = M_PI_2;                                                    
         #endif                                                                       
                                                                                      
      #elif DIM == 3                                                                  
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = grid.X2[j];                                                
         local_grid.x[3] = grid.X3[k];                                                
                                                                                      
      #endif 

      #if PHYSICS == RHD
         Get_Metric_Components(&local_grid);
      #endif 

         r   = local_grid.x[1];
         th  = local_grid.x[2];
         phi = local_grid.x[3];

         y_rr     = local_grid.gamma_con[0][0];
         y_thth   = local_grid.gamma_con[1][1];
         y_phiphi = local_grid.gamma_con[2][2];
         y_rth    = local_grid.gamma_con[0][1];
         y_rphi   = local_grid.gamma_con[0][2];
         y_thphi  = local_grid.gamma_con[1][2];

         beta_r   = local_grid.beta_con[0];
         beta_th  = local_grid.beta_con[1];
         beta_phi = local_grid.beta_con[2];

         lapse = local_grid.lapse;
         dety  = local_grid.dety;

         rho  = B(RHO,i,j);
         pre  = B(PRE,i,j);
         vr   = B(VX1,i,j);
         vth  = B(VX2,i,j);
         vphi = B(VX3,i,j);

         Vr   = y_rr*vr + y_rth*vth + y_rphi*vphi;
         Vth  = y_rth*vr + y_thth*vth + y_thphi*vphi;
         Vphi = y_rphi*vr + y_thphi*vth + y_phiphi*vphi;
         
         VV = Vr*vr + Vth*vth + Vphi*vphi;
         W  = 1.0/sqrt(1.0 - VV);

         Ur = W*(Vr - beta_r/lapse);

         if(j == gc+1 || j == Nx2-gc-1)
         {
            Mdot = Mdot - 4.0*M_PI*rho*Ur*lapse*dety*(grid.X2p[j] - grid.X2m[j])/(4.0*M_PI*density_inf);
         }
         else
         {
            Mdot = Mdot - 4.0*M_PI*rho*Ur*lapse*dety*(grid.X2p[j] - grid.X2m[j])/(4.0*M_PI*density_inf);
         }
      }
      Mean[i] = Mdot;
      Mdot    = 0.0;
      Vol     = 0.0;
   }

   // Clean the Mdot/Mdot_Mean_Time.dat file if t = 0
   file_time   = fopen(name_time,"a");
   if(time == 0.0)
   {
      char dum1[50] = "cat /dev/null > ";
      strcat(dum1," ");
      strcat(dum1,name_time);
      system(dum1);
      Mdot_Mean = 0.0;
   }

   // Computes both Mdot_Mean all over the domain
   // and Mdot at the outer horizon
   for(int i = gc + 2; i <= Nx1-gc-2; i++)
   {
      local_grid.x[1] = grid.X1[i];                                                
      local_grid.x[2] = 0.0;                                                       
      local_grid.x[3] = 0.0;                                                       
      #if COORDINATES == SPHERICAL                                                 
      local_grid.x[2] = M_PI_2;                                                    
      #endif                                                                       

      #if PHYSICS == RHD
      Get_Metric_Components(&local_grid);
      #endif 

      r   = local_grid.x[1];
      th  = local_grid.x[2];
      phi = local_grid.x[3];

      y_rr     = local_grid.gamma_con[0][0];
      y_thth   = local_grid.gamma_con[1][1];
      y_phiphi = local_grid.gamma_con[2][2];
      y_rth    = local_grid.gamma_con[0][1];
      y_rphi   = local_grid.gamma_con[0][2];
      y_thphi  = local_grid.gamma_con[1][2];

      beta_r   = local_grid.beta_con[0];
      beta_th  = local_grid.beta_con[1];
      beta_phi = local_grid.beta_con[2];

      lapse = local_grid.lapse;
      dety  = local_grid.dety;
                                                                                      
      Mdot += Mean[i];

      if(local_grid.x[1] > 1.0 + sqrt(1.0 - pow(Black_Hole_Spin,2.0)) && local_grid.x[1] <     1.0 + sqrt(1.0 - pow(Black_Hole_Spin,2.0)) + (grid.X1p[i] - grid.X1m[i]))
      {                                                                         
         Mdot_Acc = Mean[i];                                                    
      }   
      else
      {
         Mdot_Acc = Mean[gc+2];
      }
   }

   // Computes the difference between Mdot_mean(n) and Mdot_mean(n-1)
   Vol   = Nx1-2*gc-4;
   Delta = Mdot/Vol - Mdot_0;

   if(Delta > 0 && plus == FALSE)
   {
      Mdot_Max = Mdot/Vol;
      plus     = TRUE;
      minus    = FALSE;
      count++;
   }

   if(Delta < 0 && minus == FALSE)
   {
      Mdot_Min = Mdot/Vol;
      plus     = FALSE;
      minus    = TRUE;
      count++;
   }

   #if MDOT_END == TRUE
   if(time == 0.0)
   {
      restart_file = FALSE;
   }

   if(restart_simulation == FALSE && fabs(1.0 - Mdot_Min/Mdot_Max) < 1.0e-02 && restart_file == FALSE)
   {
      tprint       = time;
      restart_file = TRUE;
      printf("\n");
      printf("\n");
      printf("Mass accretion rate oscillation error (ΔM/M) drops below %e. Restart file created.\n",1.0e-02);
      printf("Simulation continues until ΔM/M < %e.\n",MDOT_ERR);
      printf("\n");
   }
   else if(Black_Hole_Spin >= 0.0 && fabs(1.0 - Mdot_Min/Mdot_Max) < MDOT_ERR && count > 10)
   {
      tmax   = time;
      tprint = tmax;
      Mdot_end  = TRUE;
      printf("\n");
      printf("\n");
      printf("Error below %e. Termination. Count = %d\n",MDOT_ERR,count);
      printf("\n");
   }
   #else
   Mdot_end = FALSE;
   #endif

   Mdot_Mean_Time = (Mdot_Max + Mdot_Min)/2.0;
   if(time == 0.0)
   {
      Mdot_Mean_Time = Mdot/Vol;
   }

   if(time >= Mdot_tprint)
   {
      fprintf(file_time,"%e %e %e %e %e\n",time,Mdot/Vol,Mdot_Acc,Mdot_Mean_Time,fabs(1.0 - Mdot_Min/Mdot_Max));
      Mdot_tprint = Mdot_tprint + MDOT_TIME;
   }
   fclose(file_time);

   if(time >= tprint)
   {
      file_domain = fopen(name_domain,"w");
      fprintf(file_domain,"###############PARAM###############\n");
      fprintf(file_domain,"%e \n",time);
      fprintf(file_domain,"%d \n",Nx1-2*gc+1);
      fprintf(file_domain,"###################################\n");

      for(int i = gc; i <= Nx1-gc; i++)
      {
         fprintf(file_domain,"%e %e\n",grid.X1[i],Mean[i]);
      }

      printf("Computed the mass accretion rate: %s\n",name_domain); 
      fclose(file_domain);
   }

   Mdot_0 = Mdot/Vol;
}
