/*                                                                              
 * File Name : user_param.h                                                     
 * Description : aztekas user parameters header file for Shock-Tube             
 * Creation Date : 27-09-2019                                                   
 * Last Modified : 28-10-2019 15:51:04
 * Created By : Alejandro Aguayo-Ortiz                                          
 */ 

#include"macros.h"

#define POLAR              TRUE

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define periodic_x2        TRUE

#define RECONST            MC
#define FLUX               HLLC
#define GRID               UNIFORM
#define PRINT_EVOLV        TRUE

double density_0, pressure_0, r_dot_0, phi_dot_0;
