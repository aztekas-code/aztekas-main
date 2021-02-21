/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Relativistic Kelvin-Helmholtz
 * Creation Date : 27-09-2019
 * Last Modified : 28-09-2019 09:47:25
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define outflow_x2max      TRUE
#define outflow_x2min      TRUE

#define RECONST            MINMOD
#define FLUX               HLL
#define GRID               UNIFORM
#define PRINT_EVOLV        TRUE

double rhotl, ptl, vx1tl, vx2tl;
double rhotr, ptr, vx1tr, vx2tr;
double rhobl, pbl, vx1bl, vx2bl;
double rhobr, pbr, vx1br, vx2br;

double x_0;
