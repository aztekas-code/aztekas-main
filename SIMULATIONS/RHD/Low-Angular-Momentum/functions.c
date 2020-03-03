/* 
 *  aztekas functions module
 *  Date of creation: 27-11-2019 10:44:40
 *  author: Alejandro Aguayo-Ortiz 
 */
#include"main.h"

double ftheta(double th)
{
   return 1.0 - fabs(cos(th));
}

double gtheta(double th)
{
   return 1.0 - delta*(1.0 - pow(sin(th),2.0));
}
