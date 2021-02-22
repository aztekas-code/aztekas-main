# Tolman-Oppenheimer-Volkoff Local Files - AZTEKAS

This the local file repository for the Tolman-Oppenheimer-Volkoff (TOV)
module for the AZTEKAS code.

The TOV module needs to be specified as:
   -> TOV=True
in the local Makefile.

This variable will set to "TRUE" the Ordinary Differential Equation (ODE)
system solver inside AZTEKAS, by enabling the ODE_Integration() function
(/aztekas-main/Src/integration/equation-solver.c). This function is defined
in the /aztekas-main/Src/integration/ODE/integration.c file, where it is already
coded the "test_module()" function.

In order to test that the TOV module works correctly, the user must create a C
file, and define the

#include"main.h"

int test_module()
{
}

NOTE: Remember to define the function in the user_param.h

Inside that function you can either print the stored values of vector U:
   -> Print_Values()
   -> Print_Values("file_id") or Print_Values()
which will be stored in the corresponding "outputdirectory" or print a
"Hello world!" to the terminal.

If either the file or the "Hello world!" are printed, that means that the TOV
module works just fine.

You can use this function for the implementation of TOV system solving, create
another function. If you do this, you must remember to call that function inside
the ODE_Integration().

# TIPS FOR IMPLEMENTING TOV solver

The solution vector is U, and its called using the following notation:
   - U(n,i), U(n,i,j) and U(n,i,j,k),
depending if you are workin with DIM = {1,2,3}, respectively, and where
where
   -> n runs from 0 to eq, where eq is the number of equations 
                          (for RHD and DIM =1, eq = 3)
   -> i,j,k runs from 0 to Nx1,Nx2,Nx3, respectively, which are defined in
     the parameter file (tov.param).

If you want to access the value of the radius (first spherical coordinate,
define COORD=Spherical in Makefile), use the variable:
   -> grid.X1[i].

The steps to follow:
   
   1. In the initial.c local file, you can define the value of U at r = 0.
   2. Then, within the test_module() function (or a user defined function)
      perform the integration forward in the domain. 
   3. You should compute the right-hand-side of you system of equations:
         dU/dr = F(U,r),
      which can, in principle, depends on U and r.
   4. Finally, integrate using Runge_Kutta(&rk,order), where rk is an already
      defined structure that contains the variables:
         rk.u0 - Zero step
         rk.u1 - First step
         rk.u2 - Second step
         rk.u3 - Third step
         rk.u4 - Fourth step
         rk.h  - Step size
         rk.f  - Right-hand-side of the system.
      In order to use this structure, you will have to declare at the beginning
      of your function as:
         -> rk_ rk;

Here is a pseudocode structure on how to implemente a 2nd Order Runge-Kutta:
http://metodos-numericos-para-ecuaciones.blogspot.com/p/metodo-de-runge-kutta-de-segundo-orden.html

   for(i)
   {
      Function_for_Computing_F(rk.f,rk.u0,r);
      Runge_Kutta(&rk,1); // Para el primer paso, se obtiene rk.u1
      Funcion_for_Computing_F(rk.f,rk.u1,r+rk.h);
      Runge_Kutta(&rk,2); // Para el segundo y final paso, se obtien rk.u2
      U(n,i) = rk.u2;
   }

You will have to specify all the variables in rk structure. May be once (like
rk.h) but, mostly, at each step i.

In AZTEKAS, you can find the implemented Runge-Kutta orders at /aztekas-main/Src/integration/runge-kutta.c (up to 3rd order).

All the files and functions created in here, could be eventually migrated to the
/aztekas-main/Src/physics/RHD/TOV/ repository.
