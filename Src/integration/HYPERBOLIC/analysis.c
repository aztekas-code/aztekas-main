/**
 * @file analysis.c
 *
 * @brief 
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 26-02-2021 - 02:30:01
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 23-02-2021 - 19:20:20
 */

#include"main.h"

void Analysis()
{
#if MDOT == TRUE
   Mass_Accretion_Rate(U);
#endif

   /*
    * Forcing termination
    */
   Alternative_Termination();
}
