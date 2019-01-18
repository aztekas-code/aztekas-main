#include<stdio.h>
#include<math.h>
#include"../Headers/matrix.h"
#include"../Headers/main.h"
    
int funct_A(double *a, double *uu)
{
   int i, j;
   double R;
   double n, p, u=0, v=0, w=0;
   n = uu[0];
   p = uu[1];
   u = uu[2];
   R = sqrt(x1*x1 + x2*x2);
   if(dim >= 2){v = uu[3];}
   if(dim == 3){w = uu[4];}
    
   for(i = 0; i < eq; i++)
   {
      for(j = 0; j < eq; j++)
      {
         if(i == 0 && j == 0)
         {
            a[i*eq + j] = 1;
         }
         else if(i == 0 && j == 1)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 0 && j == 2)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 0 && j == 3)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 0 && j == 4)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 1 && j == 0)
         {
            a[i*eq + j] = ((K-1)*pow(w,2.0)+(K-1)*pow(v,2.0)+(K-1)*pow(u,2.0))/2;
         }
         else if(i == 1 && j == 1)
         {
            a[i*eq + j] = K-1;
         }
         else if(i == 1 && j == 2)
         {
            a[i*eq + j] = (1-K)*u;
         }
         else if(i == 1 && j == 3)
         {
            a[i*eq + j] = (1-K)*v;
         }
         else if(i == 1 && j == 4)
         {
            a[i*eq + j] = (1-K)*w;
         }
         else if(i == 2 && j == 0)
         {
            a[i*eq + j] = -u/n;
         }
         else if(i == 2 && j == 1)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 2 && j == 2)
         {
            a[i*eq + j] = 1/n;
         }
         else if(i == 2 && j == 3)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 2 && j == 4)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 3 && j == 0)
         {
            a[i*eq + j] = -v/n;
         }
         else if(i == 3 && j == 1)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 3 && j == 2)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 3 && j == 3)
         {
            a[i*eq + j] = 1/n;
         }
         else if(i == 3 && j == 4)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 4 && j == 0)
         {
            a[i*eq + j] = -w/n;
         }
         else if(i == 4 && j == 1)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 4 && j == 2)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 4 && j == 3)
         {
            a[i*eq + j] = 0;
         }
         else if(i == 4 && j == 4)
         {
            a[i*eq + j] = 1/n;
         }
      }
   }
     
   return 0;
}
