##############################################
########### Makefile for Stiff EoS ###########
##############################################
#---------------------------------------------
# You can add here other functions helpful for 
# the equation of state manipulation
#---------------------------------------------

$(shell gfortran -c $(AZTPATH)/Src/physics/EOS/$(EOS)/helmholtz.f90)
FORTRAN = helmholtz.o

PHYSICS += $(AZTPATH)/Src/physics/EOS/$(EOS)/eos.c
