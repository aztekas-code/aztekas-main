##################################################
#####          MAKEFILE FOR AZTEKAS          #####
##################################################
#-------------------------------------------------
# From the list below, please activate/deactivate 
# the options that apply to your run. If you 
# modify any of these options, make sure that you 
# recompile the whole code by typing 
# $ make clean && make 
#-------------------------------------------------

##################################################
# PATH TO AZTEKAS
##################################################

AZTPATH = $(AZTEKAS_PATH)

##################################################
# SET NUMBER OF THREADS USED BY OMP
##################################################
# You can find the number of physical cores 
# using
# $ grep '^core id' /proc/cpuinfo |sort -u|wc -l
##################################################

OMP_NUM = 1

##################################################
# Physics (HD)
#    - HYDRO (Euler inviscid equations)
#    - Lane-Endem?
# Physics (RHD)
#    - HYDRO (Euler inviscid equations)
#    - TOV (Tolman-Oppenheimer-Volkoff, RHD only)
##################################################

PHY = RHD

HYDRO = True

##################################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
##################################################

DIM = 2

##################################################
# Metric (User, Minkowski, Schwarzschild, 
# Eddington-Finkelstein, Boyer-Lindquist or
# Kerr-Schild)
# ONLY IF PHY == RHD or RMHD
##################################################

METRIC = Minkowski

##################################################
# Equation of State (User, Ideal, Dust, Stiff)
##################################################

EOS = Ideal

##################################################
# Coordinates (Cartesian, Cylindrical, Spherical)
##################################################

COORD = Cylindrical

##################################################
# User Source Terms (True,False)
##################################################

USER_SOURCE = False

##################################################
# Integration method
##################################################

INT = Standard

##################################################
# All user *.c are written here
##################################################

USR = initial.c \
      user_boundaries.c \
		user_input.c 

##################################################
# Include a global Makefile in the /Src/ 
# directory, which has some definitions
# of parameters, inclusions of other makefiles
# and compilation orders. Change this Makefile,
# and recursive Makefiles only if you know what
# your doing.
##################################################
include $(AZTPATH)/Src/Makefile
