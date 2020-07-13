/**
 * @file /RHD/surface_volume.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 */

#include"main.h"

void Surface_Volume()
{
   int i, j, k;
   double volimjpkp, volijpkp, volipjpkp;
   double volimjkp,  volijkp,  volipjkp;
   double volimjmkp, volijmkp, volipjmkp;
   double volimjpk,  volijpk,  volipjpk;
   double volimjk,   volijk,   volipjk;
   double volimjmk,  volijmk,  volipjmk;
   double volimjpkm, volijpkm, volipjpkm;
   double volimjkm,  volijkm,  volipjkm;
   double volimjmkm, volijmkm, volipjmkm;
   double dV;
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
         // ip jp k
         vol.x[1]    = grid.X1p[i];
         vol.x[2]    = grid.X2p[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volipjpk = vol.dety;
 
         // i jp k
         vol.x[1]    = grid.X1[i];
         vol.x[2]    = grid.X2p[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volijpk = vol.dety;

         // im jp k
         vol.x[1]    = grid.X1m[i];
         vol.x[2]    = grid.X2p[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volimjpk = vol.dety;
 
         // ip j k
         vol.x[1]    = grid.X1p[i];
         vol.x[2]    = grid.X2[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volipjk = vol.dety;
 
         // i j k
         vol.x[1]    = grid.X1[i];
         vol.x[2]    = grid.X2[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volijk = vol.dety;

         // im j k
         vol.x[1]    = grid.X1m[i];
         vol.x[2]    = grid.X2[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volimjk = vol.dety;
 
         // ip jm k
         vol.x[1]    = grid.X1p[i];
         vol.x[2]    = grid.X2m[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volipjmk = vol.dety;
 
         // i jm k
         vol.x[1]    = grid.X1[i];
         vol.x[2]    = grid.X2m[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volijmk = vol.dety;

         // im jm k
         vol.x[1]    = grid.X1m[i];
         vol.x[2]    = grid.X2m[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         #endif
         vol.x[3]    = 0.0;

         Get_Metric_Components(&vol);
         volimjmk = vol.dety;

         //dV = (      (volimjpk + 4.0*volijpk + volipjpk)/6.0 \
               + 4.0*( volimjk + 4.0* volijk +  volipjk)/6.0 \
               +     (volimjmk + 4.0*volijmk + volipjmk)/6.0)/6.0;
         //dV = (  4.0*volimjpk +  8.0*volijpk + 4.0*volipjpk \
               + 8.0*volimjk  + 16.0*volijk  + 8.0*volipjk  \
               + 4.0*volimjmk +  8.0*volijmk + 4.0*volipjmk)/64;
         dV = (  1.0*volimjpk +  4.0*volijpk + 1.0*volipjpk \
               + 4.0*volimjk  + 16.0*volijk  + 4.0*volipjk  \
               + 1.0*volimjmk +  4.0*volijmk + 1.0*volipjmk)/36.0;
         //dV = volijk;

         vol.x[1]    = grid.X1[i];
         surf_p.x[1] = grid.X1p[i];
         surf_m.x[1] = grid.X1m[i];
         vol.x[2]    = grid.X2[j];
         surf_p.x[2] = grid.X2[j];
         surf_m.x[2] = grid.X2[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         surf_p.x[2] = M_PI_2;
         surf_m.x[2] = M_PI_2;
         #endif
         vol.x[3]    = 0.0;
         surf_p.x[3] = 0.0;
         surf_m.x[3] = 0.0;
 
         Get_Metric_Components(&vol);
         Get_Metric_Components(&surf_p);
         Get_Metric_Components(&surf_m);
 
         //S1p(i,j) = surf_p.dety/((surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0);
         //S1m(i,j) = surf_m.dety/((surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0);
         S1p(i,j) = surf_p.dety/dV;
         S1m(i,j) = surf_m.dety/dV;
         //S1p(i,j) = surf_p.dety/vol.dety;
         //S1m(i,j) = surf_m.dety/vol.dety;

         //printf("\n");
         //printf("%e %e %e\n",dV,(surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0,volijk);
         //printf("\n");
 
         vol.x[1]    = grid.X1[i];
         surf_p.x[1] = grid.X1[i];
         surf_m.x[1] = grid.X1[i];
         vol.x[2]    = grid.X2[j];
         surf_p.x[2] = grid.X2p[j];
         surf_m.x[2] = grid.X2m[j];
         #if POLAR == TRUE
         vol.x[2]    = M_PI_2;
         surf_p.x[2] = M_PI_2;
         surf_m.x[2] = M_PI_2;
         #endif
         vol.x[3]    = 0.0;
         surf_p.x[3] = 0.0;
         surf_m.x[3] = 0.0;
 
         Get_Metric_Components(&vol);
         Get_Metric_Components(&surf_p);
         Get_Metric_Components(&surf_m);
 
         //S2p(i,j) = surf_p.dety/((surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0);
         //S2m(i,j) = surf_m.dety/((surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0);
         S2p(i,j) = surf_p.dety/dV;
         S2m(i,j) = surf_m.dety/dV;
         //S2p(i,j) = surf_p.dety/vol.dety;
         //S2m(i,j) = surf_m.dety/vol.dety;

         //printf("\n");
         //printf("%e %e %e\n",dV,(surf_p.dety + 4.0*vol.dety + surf_m.dety)/6.0,volijk);
         //printf("\n");
 
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
