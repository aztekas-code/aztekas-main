#include"main.h"

void Surface_Volume()
{
   int i, j, k;
   gauge_ vol, surf_p, surf_m;

#if DIM == 1

   for(i = 0; i <= Nx1; i++)
   {
      vol.x[1]    = grid.X1[i];
      surf_p.x[1] = grid.X1p[i];
      surf_m.x[1] = grid.X1m[i];
      vol.x[2]    = M_PI_2;
      surf_p.x[2] = M_PI_2;
      surf_m.x[2] = M_PI_2;
      vol.x[3]    = 0.0;
      surf_p.x[3] = 0.0;
      surf_m.x[3] = 0.0;

      Get_Metric_Components(&vol);
      Get_Metric_Components(&surf_p);
      Get_Metric_Components(&surf_m);

      S1p(i) = surf_p.dety/vol.dety;
      S1m(i) = surf_m.dety/vol.dety;
   }
   
#elif DIM == 2 || DIM == 4

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         vol.x[1]    = grid.X1[i];
         surf_p.x[1] = grid.X1p[i];
         surf_m.x[1] = grid.X1m[i];
         vol.x[2]    = grid.X2[j];
         surf_p.x[2] = grid.X2[j];
         surf_m.x[2] = grid.X2[j];
         vol.x[3]    = 0.0;
         surf_p.x[3] = 0.0;
         surf_m.x[3] = 0.0;
 
         Get_Metric_Components(&vol);
         Get_Metric_Components(&surf_p);
         Get_Metric_Components(&surf_m);
 
         S1p(i,j) = surf_p.dety/vol.dety;
         S1m(i,j) = surf_m.dety/vol.dety;
 
         vol.x[1]    = grid.X1[i];
         surf_p.x[1] = grid.X1[i];
         surf_m.x[1] = grid.X1[i];
         vol.x[2]    = grid.X2[j];
         surf_p.x[2] = grid.X2p[j];
         surf_m.x[2] = grid.X2m[j];
 
         Get_Metric_Components(&vol);
         Get_Metric_Components(&surf_p);
         Get_Metric_Components(&surf_m);
 
         S2p(i,j) = surf_p.dety/vol.dety;
         S2m(i,j) = surf_m.dety/vol.dety;
      }
   }

#elif DIM == 3

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         for(k = 0; k <= Nx3; k++)
         {
            S1p(i,j,k) = 1.0;
            S1m(i,j,k) = 1.0;
            S2p(i,j,k) = 1.0;
            S2m(i,j,k) = 1.0;
            S3p(i,j,k) = 1.0;
            S3m(i,j,k) = 1.0;
         }
      }
   }

#endif
}
