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

OMP_NUM = 3

###############################################
# Physics (hd, mhd, rhd, rmhd)
###############################################
# You can find the number of threads using
# $ lscpu
# or
# $ nproc
# Be careful of using only the number of 
# PHYSICAL cores
###############################################

###############################################
# Physics (HD, MHD, RHD, RMHD)
###############################################

PHY = HD

###############################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 2

########################################################################################
# Metric (Minkowski, Schwarzschild, Eddington-Finkelstein, Boyer-Lindquist, Kerr-Schild)
########################################################################################

## ONLY IF PHY == rhd or rmhd

METRIC = Minkowski

###############################################
# Equation of State (ideal, dust, stiff)
###############################################

EOS = ideal

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = cartesian

###############################################
# User Source Terms (true,false)
###############################################

USER_SOURCE = true

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
