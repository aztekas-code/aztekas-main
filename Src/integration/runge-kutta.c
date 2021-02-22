/**
 * @file /integration/runge-kutta.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Runge-Kutta steps.
 * 
 * This is quite different from the standar Runge-Kutta method.
 * In the standard form, each step is called k_i, and corresponds to
 * k_i = h*F(\sum_j \propto k_{j < i}), 
 * and finally,
 * U_{n+1} = U_n + \sum_i \propto k_i.
 *
 * But in this way, we obtain different steps for U itself, without the
 * intermediate k_i. 
 */

#include"main.h"

void Runge_Kutta(rk_ *rk, int step)
{
   if(step == 1)
   {
      rk->u1 = rk->u0 + rk->h*rk->f;
   }

   if(step == 2 && rk_order == 2)
   {
      rk->u2 = 0.5*(rk->u0 + rk->u1 + rk->h*rk->f);
   }

   if(step == 2 && rk_order == 3)
   {
      rk->u2 = 0.25*(3.0*rk->u0 + rk->u1 + rk->h*rk->f);
   }

   if(step == 3)
   {
      rk->u3 = (1.0/3.0)*(3.0*rk->u0 + 2.0*rk->u2 + 2.0*rk->h*rk->f);
   }
}
