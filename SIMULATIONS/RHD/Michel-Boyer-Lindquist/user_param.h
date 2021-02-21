/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Relativistic Shock-Tube
 * Creation Date : 28-09-2019
 * Last Modified : 17-12-2019 13:10:55
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define POLAR              FALSE

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#if POLAR == FALSE
   #define reflective_x2max   TRUE
   #define reflective_x2min   TRUE
#elif POLAR == TRUE
   #define periodic_x2        TRUE
#endif

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM
#define PRINT_EVOLV        TRUE

#define Black_Hole_Mass    1.0
#define Black_Hole_Spin    0.0

double density_0, pressure_0, velocity_0;
