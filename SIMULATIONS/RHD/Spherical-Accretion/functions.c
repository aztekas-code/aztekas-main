/**
 * @file functions.c
 * 
 * @brief Useful functions
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 24-10-2020 - 00:29:37
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 28-04-2020
 */

#include"main.h"

void Mass_Accretion_Rate(double *B)
{
   /* 
    * File variables
    */
   FILE *file_time, *file_acc,\
        *file_domain, *file_theta;
   char name_time[100], name_acc[100],\
        name_domain[100], name_theta[100];
   char ext[50], dat[50];
   int num;

   /*
    * Local grid
    */
   gauge_ local_grid;

   /*
    * Important local variables
    */
   #if COORDINATES == SPHERICAL
   double time, r, th, thm, phi;
   double dT, dr, dth, dphi;
   double y_rr, y_thth, y_phiphi;
   double y_rth, y_rphi, y_thphi;
   double beta_r, beta_th, beta_phi;
   double lapse, dety;
   double rho2, r_plus;
   double rho, pre, vr, vth, vphi;
   double Vr, Vth, Vphi, VV, W;
   double UT, Ur, Ut, Uphi;
   double UhT, Uhr, Uht, Uhphi;
   double UhatT, Uhatr, Uhatt;
   double Uhracc, Uhtacc, racc;
   double Urmean[Nx1+1], Utmean[Nx1+1], dV[Nx1+1];
   double Uhtr1[Nx2+1], Uhtr2[Nx2+1];
   double Uhtr3[Nx2+1], Uhtr4[Nx2+1], Uhtr5[Nx2+1];
   double Mdot, Min, Mej, Vol, dVol;
   double Mdotmean[Nx1+1];
   double Mdot_Acc = 0.0;
   double Mdot_Mean_Time = 0.0;
   double Delta;
   #endif

   /*
    * Initialize time and event horizon r_+
    */
   local_grid.x[0] = grid.time;
   time            = local_grid.x[0];
   r_plus          = 1.0 + sqrt(1.0 - pow(Black_Hole_Spin,2.0));

   /*
    * Initialize average values of Mdot, Uhatr and Uhatt
    */
   Mdot  = 0.0;
   Uhatr = 0.0;
   Uhatt = 0.0;
   dVol  = 0.0;

   /*
    * Only if time greater than the Mdot_tprint time
    * or when time == tmax, do the calculations
    */
   if(time >= Mdot_tprint || time >= tmax)
   {
      /*
       * Number of the file to print and extension
       */
      num = itprint;
      strcpy(ext,".dat");
      snprintf(dat,8,"%d",num);

      /*
       * Path and name of Mdot_Mean_Time.dat
       * Mdot_Acc_Time.dat, Mdot_Mean_Domain_<num>.dat
       * and Utheta_<num>.dat
       */
      strcpy(name_time,outputdirectory);
      strcat(name_time,"Analysis/");
      strcat(name_time,"Average_vs_Time");
      strcat(name_time,ext);

      strcpy(name_acc,outputdirectory);
      strcat(name_acc,"Analysis/");
      strcat(name_acc,"Acc_vs_Time");
      strcat(name_acc,ext);

      strcpy(name_domain,outputdirectory);
      strcat(name_domain,"Analysis/");
      strcat(name_domain,"Value_vs_r_");
      strcat(name_domain,dat);
      strcat(name_domain,ext);

      strcpy(name_theta,outputdirectory);
      strcat(name_theta,"Analysis/");
      strcat(name_theta,"Utheta_");
      strcat(name_theta,dat);
      strcat(name_theta,ext);

      /*
       * Integration loop avoiding \theta = 0 and \theta_max
       */
      for(int i = gc; i <= Nx1-gc; i++)
      {
         for(int j = gc+1; j <= Nx2-gc-0; j++)
         {
            /*
             * Fill local_grid structure
             */
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

            /*
             * Get 3+1 gauge metric components
             */
         #if PHYSICS == RHD
            Get_Metric_Components(&local_grid);
         #endif 

            /*
             * Obtain r, th, phi, dr, dth, alpha, beta^i, y^ij, dety and rho2
             */
            r   = local_grid.x[1];
            th  = local_grid.x[2];
            phi = local_grid.x[3];
            dr  = grid.X1p[i] - grid.X1m[i];
            dth = grid.X2p[j] - grid.X2m[j];

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

            rho2 = r*r + Black_Hole_Spin*Black_Hole_Spin*cos(th)*cos(th);

            /*
             * Primitive variables rho, P, v_r, v_th, v_phi
             */
            rho  = B(RHO,i,j);
            pre  = B(PRE,i,j);
            vr   = B(VX1,i,j);
            vth  = B(VX2,i,j);
            vphi = B(VX3,i,j);

            /*
             * v^r, v^th, v^phi
             */
            Vr   = y_rr*vr + y_rth*vth + y_rphi*vphi;
            Vth  = y_rth*vr + y_thth*vth + y_thphi*vphi;
            Vphi = y_rphi*vr + y_thphi*vth + y_phiphi*vphi;
            
            /*
             * v_i*v^i and Lorentz factor W
             */
            VV = Vr*vr + Vth*vth + Vphi*vphi;
            W  = 1.0/sqrt(1.0 - VV);

            /*
             * Four-velicity U^\mu
             */
            UT   = W/lapse;
            Ur   = W*(Vr   - beta_r/lapse);
            Ut   = W*(Vth  - beta_th/lapse);
            Uphi = W*(Vphi - beta_phi/lapse);

            /*
             * Four-velocity in orthonormal coordinates U^\hat{\mu}
             */
            UhT   = (1.0/sqrt(1.0 + 2.0*r/rho2))*UT;  
            Uhr   = sqrt(1.0 + 2.0*r/rho2)*(Ur + (2.0*r/(rho2 + 2.0*r))*UT \
                    - Black_Hole_Spin*pow(sin(th),2.0)*Uphi);
            Uht   = sqrt(rho2)*Ut;
            Uhphi = sqrt(rho2)*sin(th)*Uphi;

            /*
             * Composite trapezoidal integration for Mdot U^\hat{r} 
             * and U^\hat{th}
             */
            if(j == gc+1 || j == Nx2-gc-0)
            {
               Mdot  = Mdot  - 4.0*M_PI*rho*Ur*lapse*dety*dth/((1.0/density_0)*M_PI)/2.0;
               Uhatr = Uhatr + 4.0*M_PI*Uhr*lapse*dety*dth/2.0;
               Uhatt = Uhatt + 4.0*M_PI*Uht*lapse*dety*dth/2.0;
               Vol   = Vol   + 4.0*M_PI*lapse*dety*dth/2.0;
            }
            else
            {
               Mdot  = Mdot  - 4.0*M_PI*rho*Ur*lapse*dety*dth/((1.0/density_0)*M_PI);
               Uhatr = Uhatr + 4.0*M_PI*Uhr*lapse*dety*dth;
               Uhatt = Uhatt + 4.0*M_PI*Uht*lapse*dety*dth;
               dVol  = dVol  + 4.0*M_PI*lapse*dety*dth;
            }

            /*
             * Store the values of U^\hat{th} for
             * r = (r_+, 2r_+, R/4, R/2, 3R/4) for each th[j],
             * with R = rmax = x1max
             */
            if(r >= r_plus && r < r_plus + dr)
            {                                                                         
               Uhtr1[j] = Uht;
            }
            else if(r >= 2.0*r_plus && r < 2.0*r_plus + dr)
            {                                                                         
               Uhtr2[j] = Uht;
            }
            else if(r >= x1max/4.0 && r < x1max/4.0 + dr)
            {                                                                         
               Uhtr3[j] = Uht;
            }
            else if(r >= x1max/2.0 && r < x1max/2.0 + dr)
            {                                                                         
               Uhtr4[j] = Uht;
            }
            else if(r >= 3.0*x1max/4.0 && r < 3.0*x1max/4.0 + dr)
            {                                                                         
               Uhtr5[j] = Uht;
            }
         }

         /*
          * Store values of Mdot Uhatr and Uhatt angular integrations
          * for each radius r[i]
          */
         Mdotmean[i] = Mdot;
         Urmean[i]   = Uhatr;
         Utmean[i]   = Uhatt;
         dV[i]       = dVol;

         /*
          * Values at the event horizon r_+ or at the fifth cell,
          * in case r_+ < rmin
          */
         if(r >= r_plus && r < r_plus + dr)
         {                                                                         
            Mdot_Acc = Mdotmean[i]; 
            Uhracc   = Urmean[i]/dV[i];
            Uhtacc   = Utmean[i]/dV[i];
            racc     = r;
         }   
         else if(r_plus < x1min)
         {
            Mdot_Acc = Mdotmean[gc+5];
            Uhracc   = Urmean[gc+5]/dV[gc+5];
            Uhtacc   = Utmean[gc+5]/dV[gc+5];
            racc     = grid.X1[gc+5];
         }

         /*
          * Re-initialize to zero the dummy variables
          */
         Mdot  = 0.0;
         Uhatr = 0.0;
         Uhatt = 0.0;
         dVol  = 0.0;
      }

      /*
       * Open and clean the Analysis/Mdot_Mean_Time.dat 
       * and Analysis/Mdot_Acc_Time.dat if file if t = 0
       */
      file_time   = fopen(name_time,"a");
      file_acc    = fopen(name_acc,"a");
      if(time == 0.0)
      {
         char dum1[50] = "cat /dev/null > ";
         strcat(dum1," ");
         strcat(dum1,name_time);
         system(dum1);
         char dum2[50] = "cat /dev/null > ";
         strcat(dum2," ");
         strcat(dum2,name_acc);
         system(dum2);
         Mdot_Mean = 0.0;
         Mdot_Max = 1.0;
         Mdot_Min = 1.0;
         fprintf(file_time,"#time,Mdot,Uhatr,Uhatth,errMdot\n");
         fprintf(file_acc,"#time,Mdot_acc,Uhatr_acc,Uhatth_acc\n");
      }

      /*
       * Computes the average values all over the domain for
       * Mdot, Uhatr and Uhatth
       */
      Vol   = Nx1-2*gc-9; // Number of radial points
      for(int i = gc + 5; i <= Nx1-gc-5; i++)
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
         dr  = grid.X1p[i] - grid.X1m[i];

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

         /*
          * Average values
          */
         Mdot  += Mdotmean[i];
         Uhatr += Urmean[i]*dr;
         Uhatt += Utmean[i]*dr;
         dVol  += dV[i]*dr;
      }

      Mdot  = Mdot/Vol;
      Uhatr = Uhatr/dVol;
      Uhatt = Uhatt/dVol;

      /************************************************************************
       * This calculations are needed in order to monitorize the oscillations
       * in the average of Mdot
       ************************************************************************/

      /*
       * Mdot_0 is a global variable in which the value of the last max or min
       * in the oscillation is stored. Delta is the difference between that
       * point and the current point.
       */
      Delta = Mdot - Mdot_0;

      /*
       * If Delta > 0, the derivative is positive. If plus is FALSE, and
       * hence, minus is TRUE, it means that the last derivative was negative
       * and so now has a change of sign. Which means we have a minimum.
       * And viceversa with Delta < 0
       */
      if(Delta > 0 && plus == FALSE)
      {
         Mdot_Min = Mdot;
         plus     = TRUE;
         minus    = FALSE;
         count++;
      }

      if(Delta < 0 && minus == FALSE)
      {
         Mdot_Max = Mdot;
         plus     = FALSE;
         minus    = TRUE;
         count++;
      }

      /*
       * If MDOT_END == TRUE ends the simulations by considering the size of the 
       * oscillations in Mdot
       */
      #if MDOT_END == TRUE
      if(time == 0.0)
      {
         restart_file = FALSE;
      }

      /*
       * If the difference between the oscillations is less than 10⁻²,
       * send a firts alert. If the differences drops below the MDOT_ERR 
       * defined in paramfile, then tmax and tprint are equal to the current
       * time. So the simulations finish.
       */
      if(restart_simulation == FALSE && fabs(1.0 - Mdot_Min/Mdot_Max) < 1.0e-02 \
         && restart_file == FALSE)
      {
         tprint       = time;
         restart_file = TRUE;
         printf("\n");
         printf("\n");
         printf("Mass accretion rate oscillation error (ΔM/M) drops below %e.",1.0e-02);
         printf("Restart file created.\n");
         printf("Simulation continues until ΔM/M < %e.\n",MDOT_ERR);
         printf("\n");
      }
      else if(Black_Hole_Spin >= 0.0 && fabs(1.0 - Mdot_Min/Mdot_Max) < MDOT_ERR \
              && count > 10)
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

      /************************************************************************
       ************************************************************************/

      /*
       * Compute the mean value of the oscillations
       */
      Mdot_Mean_Time = (Mdot_Max + Mdot_Min)/2.0;
      if(time == 0.0)
      {
         Mdot_Mean_Time = Mdot;
      }

      /*
       * Print values in Mdot_Mean_Time.dat
       */
      fprintf(file_time,"%e %e %e %e %e\n",time,\
                                           Mdot,\
                                           Uhatr,\
                                           Uhatt,\
                                           fabs(1.0 - Mdot_Min/Mdot_Max));
      Mdot_tprint = Mdot_tprint + MDOT_TIME;
      fclose(file_time);

      /*
       * Print values in Mdot_Acc_Time.dat
       */
      fprintf(file_acc,"%e %e %e %e\n",time,\
                                        Mdot_Acc,\
                                        Uhracc,\
                                        Uhtacc);
      fclose(file_acc);

      /*
       * Only when time is tprint, print the files Mdot_Mean_Domain_<num>.dat
       * and Utheta_<num>.dat
       */
      if(time >= tprint)
      {
         file_domain = fopen(name_domain,"w");
         fprintf(file_domain,"###############PARAM###############\n");
         fprintf(file_domain,"%e \n",time);
         fprintf(file_domain,"%d \n",Nx1-2*gc+1);
         fprintf(file_domain,"#r,Mdot,Uhatr,Uhatth\n");
         fprintf(file_domain,"###################################\n");

         for(int i = gc; i <= Nx1-gc; i++)
         {
            fprintf(file_domain,"%e %e %e %e\n",grid.X1[i],\
                                                Mdotmean[i],\
                                                Urmean[i]/dV[i],\
                                                Utmean[i]/dV[i]);
         }

         fclose(file_domain);

         file_theta = fopen(name_theta,"w");
         fprintf(file_theta,"###############PARAM###############\n");
         fprintf(file_theta,"%e \n",time);
         fprintf(file_theta,"%d \n",Nx2-2*gc+1);
         fprintf(file_theta,"#theta,Uhatth={r_+,2r_+,R/4,R/2,3R/4}\n");
         fprintf(file_theta,"###################################\n");

         for(int j = gc+1; j <= Nx2-gc-1; j++)
         {
            fprintf(file_theta,"%e %e %e %e %e %e\n",grid.X2[j],\
                                                     Uhtr1[j],\
                                                     Uhtr2[j],\
                                                     Uhtr3[j],\
                                                     Uhtr4[j],\
                                                     Uhtr5[j]);
         }

         fclose(file_theta);
      }

      /*
       * Store the value of the previous time
       */
      Mdot_0 = Mdot;

   }

   /*
    * Finally prints the relevant results and parameters
    */
   if(time >= tmax)
   {
      FILE *results;
      strcpy(name_time,outputdirectory);
      strcat(name_time,"../");
      strcat(name_time,"results_a");
      strcat(name_time,ext);

      results = fopen(name_time,"a");
      fprintf(results,"%e %e %e %e %e %e\n",racc,\
                                            Temp,\
                                            cs,\
                                            Black_Hole_Spin,\
                                            Mdot_Acc,
                                            Mdot_Mean_Time);
      fclose(results);

      strcpy(name_time,outputdirectory);
      strcat(name_time,"../../");
      strcat(name_time,"results");
      strcat(name_time,ext);

      results = fopen(name_time,"a");
      fprintf(results,"%e %e %e %e %e %e\n",racc,\
                                            Temp,\
                                            cs,\
                                            Black_Hole_Spin,\
                                            Mdot_Acc,\
                                            Mdot_Mean_Time);
      fclose(results);
   }
}

double rs(double K, double c_s)
{
   double n       = 1.0/(K - 1.0);
   double h_inf   = 1.0/(1.0 - n*pow(c_s,2.0));
   double Psi     = acos((3.0/(2.0*n*h_inf))*pow((n + 3.0)/(3.0*n),-1.5))/3.0;
   double h_c     = 2.0*h_inf*sqrt((n + 3.0)/(3.0*n))*sin(Psi + M_PI/6.0);
   double c_s_c2  = (h_c - 1.0)/(n*h_c);
   double v_c2    = c_s_c2/(1.0 + 3.0*c_s_c2);
   double rs      = 0.5/v_c2;

   return rs;
}
