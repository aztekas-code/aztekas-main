/**
 * @file /integration/runge-kutta.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Principal loop in the integration.
 */

#include"main.h"

void Runge_Kutta(rk_ *rk, int order)
{
   if(order == 1)
   {
      rk->u1 = rk->u0 + rk->h*rk->f;
   }

   if(order == 2)
   {
      rk->u2 = 0.5*(rk->u0 + rk->u1 + rk->h*rk->f);
   }

   if(order == 3)
   {
      rk->u3 = (1.0/3.0)*(3.0*rk->u0 + 2.0*rk->u2 + 2.0*rk->h*rk->f);
   }
}
