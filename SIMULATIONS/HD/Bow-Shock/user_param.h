/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 15:43:29
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1min      TRUE
#define outflow_x1max      TRUE
#define outflow_x2min      TRUE
#define outflow_x2max      TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               UNIFORM
#define PRINT_EVOLV        TRUE

double r_in, r_out;
