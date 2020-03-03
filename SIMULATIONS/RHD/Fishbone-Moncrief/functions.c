/* 
 *  aztekas functions module
 *  Date of creation: 27-11-2019 10:44:40
 *  author: Alejandro Aguayo-Ortiz 
 */
#include"main.h"

double rho(double x)
{
   double a      =  1.2268882902103233E+01;
   double b      =  3.0781593823388236E-02;
   double c      = -1.3596426267607620E+00;
   double Offset = -9.8034182617611407E+01;

   return exp(a + (b/x) + c*log(x)) + Offset;
}

double pre(double x)
{
   double a      =  1.0729543490704748E+01;
   double b      =  4.8257535692800071E-02;
   double c      = -2.2707279877581454E+00;
   double Offset = -1.5252713681060741E+00;

   return  exp(a + (b/x) + c*log(x)) + Offset;
}

double vx1(double x)
{
   double a = -1.0436021879029216E-01;
   double b = -5.9152587898346848E-02;
   double c =  3.0615116715089069E-02;
   double d = -3.3500918067639676E-03;

   return a + b*log(x) + c*pow(log(x),2.0) + d*pow(log(x),3.0);
}

double ftheta(double th)
{
   return 1.0 - fabs(cos(th));
}

double gtheta(double th)
{
   return 1.0 - delta*(1.0 - pow(sin(th),2.0));
}
