/**
 * @file /integration/integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

#include"main.h"

void Hyperbolic_Integration()
{
   int rk_order = 2;

   /**
    * Convert all primitive vector U to conservative Q0
    */
   Prim2Cons_All(Q0,U);
   
   for(int order = 1; order <= rk_order; order++)
   {
      Primitive_Reconstruction();
      Method_of_Lines(order);
      Cons2Prim(U,Q);
      Boundaries(U);
      U0 = U;
   }

   /**
    * Increase time by dt
    */
   grid.time = grid.time + dt;
}

/**
 * En esta función puedes realizar todas las pruebas necesarias para la
 * integración de las ecuaciones de TOV. Define todas las variables y funciones
 * que necesites en este archivo, y luego las moveremos de lugar.
 *
 * El vector solución es U, y se llama usando la notación:
 *    - U(n,i), U(n,i,j) y U(n,i,j,k), dependiendo de si está definido en
 *      un dominio de 1, 2 o 3 dimensiones, respectivamente.
 *
 * En el caso de TOV se usará DIM = 1, por lo que tenemos U(n,i).
 *
 * La variable n debe correr de 0 a eq, donde eq es el número de ecuaciones.
 * La variable i debe correr de 0 a Nx1, donde Nx1 es el número de puntos del
 * dominio. Este es un valor fijo.
 *
 * Para acceder al valor del radio r, usa la variable grid.X1[i], donde la i, 
 * nuevamente corre de 0 a Nx1, con grid.X1[0] = (r=0) y grid.X1[Nx1] =
 * (r=R_max).
 *
 * Los pasos a seguir deben ser los siguientes:
 *
 *    1. En el archivo initial.c de tu carpeta local, puedes definir el valor en 
 *    r = 0, para poder empezar a integrar. Es decir, defines U(n,0) = U_0(n)
 *
 *    2. Posteriormente, en esta función haces la integración hacia adelante en 
 *    el dominio. En principio los siguientes pasos deben de ir dentro de un 
 *    loop. Ya sea un for(i = 0; i <= Nx1; i++) o un while(r < R_max).
 *
 *    3. Debes calcular la función del lado derecho de tu sistema de ecuaciones
 *                            dU/dr = F(U,r),
 *       la cual puede depender tanto de U, como de r.
 *
 *    4. Posteriormente integras usando la función Runge_Kutta(&rk,order), donde
 *    rk es una estructura que contiene las variables
 *          rk.u0 - Paso cero
 *          rk.u1 - Primer paso
 *          rk.u2 - Segundo paso
 *          rk.u3 - Tercer paso
 *          rk.u4 - Cuarto paso
 *          rk.h  - El intervalo de paso entre punto y punto del dominio.
 *          rk.f  - El lado derecho de la ecuación.
 *    Esta estructura se debe declarar al inicio de la función como rk_ rk;
 *
 *    Este es un pseudocódigo de ejemplo para Runge-Kutta orden 2:
 *    http://metodos-numericos-para-ecuaciones.blogspot.com/p/metodo-de-runge-kutta-de-segundo-orden.html
 *
 *    for(i)
 *    {
 *       Funcion_Que_Calcula_F(rk.f,rk.u0,r);
 *       Runge_Kutta(&rk,1); // Para el primer paso, se obtiene rk.u1       
 *       Funcion_Que_Calcula_F(rk.f,rk.u1,r+rk.h);
 *       Runge_Kutta(&rk,2); // Para el segundo y final paso, se obtien rk.u2
 *       U(n,i) = rk.u2;
 *    }
 *
 *    Para orden 3, que es lo que hay implementado seguir la siguiente página.
 *    http://metodos-numericos-para-ecuaciones.blogspot.com/p/metodo-de-runge-kutta-de-tercer-orden.html
 */
void ODE_Integration()
{
}
