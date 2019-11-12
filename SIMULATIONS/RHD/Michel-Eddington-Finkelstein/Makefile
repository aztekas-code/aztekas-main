######################################
######## MAKEFILE FOR AZTEKAS ########
######################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------

###############################################
# PATH TO AZTEKAS
###############################################

AZTPATH = $(AZTEKAS_PATH)

###############################################
# SET NUMBER OF THREADS USED BY OMP
###############################################
# You can find the number of physical cores 
# using
# $ grep '^core id' /proc/cpuinfo |sort -u|wc -l
###############################################

OMP_NUM = 1

###############################################
# Physics (HD, MHD, RHD, RMHD)
###############################################

PHY = RHD

###############################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 1

########################################################################################
# Metric (Minkowski, Schwarzschild, Eddington-Finkelstein, Boyer-Lindquist, Kerr-Schild)
########################################################################################

## ONLY IF PHY == rhd or rmhd

METRIC = Eddington-Finkelstein

###############################################
# Equation of State (User, Ideal, Dust, Stiff)
###############################################

EOS = Ideal

###############################################
# Coordinates (Cartesian, Cylindrical, Spherical)
###############################################

COORD = Spherical

###############################################
# User Source Terms (true,false)
###############################################

USER_SOURCE = false

###############################################
# Integration method
###############################################

INT = standard

###############################################

###############################################
# All user *.c are written here
###############################################

USR = initial.c \
		user_boundaries.c \
		user_input.c 

###############################################
# Include a global Makefile in the /Src/ 
# directory, which include some definitions
# of parameters, inclusions of other makefiles
# and compilation orders
###############################################
include $(AZTPATH)/Src/Makefile
