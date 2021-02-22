# Tolman-Oppenheimer-Volkoff Local Files - AZTEKAS

This the local file repository for the Tolman-Oppenheimer-Volkoff (TOV)
module for the AZTEKAS code.

The TOV module needs to be specified as:
   TOV=True
in the local Makefile.

This variable will set to "TRUE" the Ordinary Differential Equation (ODE)
system solver inside AZTEKAS, by enabling the ODE_Integration() function
(/aztekas-main/Src/integration/equation-solver.c). This function is defined
in the /aztekas-main/Src/integration/integration.c file, where it is already
coded the "hello_world_tov()".

In order to test that the TOV module works correctly.
