# Low-Angular Momentum Test - AZTEKAS

In this repository you can find the local AZTEKAS files for the Low-Angular
Momentum accretion problem. Specifically the Model 1b of Mach et al. (2008) 
https://ui.adsabs.harvard.edu/abs/2018CQGra..35i5005M/abstract

Parameters:

x1max = 300 # Define at user_input.c
x1min = 1.5 # Define at user_input.c
l_0   = 6.657 # Define at user_param.h
c_s   = 0.070710 # Computed in user_input.c
Theta = 3.0226e-03 # Define at sim
r_s   = 11.867224 # Computed in user_input.c

Run as:

./sim 1.6666666666
