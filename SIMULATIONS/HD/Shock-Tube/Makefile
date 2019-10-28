<<<<<<< HEAD
#####################################
#### Define dimension parameters ####
#####################################
=======
######################################
######## MAKEFILE FOR AZTEKAS ########
######################################
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
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
<<<<<<< HEAD
###############################################

PARAM += -DOMP_NUM=3

###############################################
# Physics (hd, mhd, rhd, rmhd)
###############################################

PHY = hd
=======
# You can find the number of threads using
# $ lscpu
# or
# $ nproc
# Be careful of using only the number of 
# PHYSICAL cores
###############################################

OMP_NUM = 1

###############################################
# Physics (HD, MHD, RHD, RMHD)
###############################################

PHY = HD
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

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
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 1

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = cartesian

###############################################
# User Source Terms (true,false)
###############################################

USER_SOURCE = false

###############################################
# Integration method
###############################################

INT = standard

###############################################

<<<<<<< HEAD
=======
###############################################
# All user *.c are written here
###############################################

>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
USR = boundaries.c \
		initial.c \
		user_input.c 

<<<<<<< HEAD
include $(AZTPATH)/Src/Makefile

###############################################
# Compilation
###############################################

FLAGS = -Ofast -lm -fopenmp

COMPILER = gcc
AZT_HEAD = $(AZTPATH)/Src/include
HEADERS = -I$(AZT_HEAD) -I.
EXEC = aztekas

all: 
	@echo ""
	@echo "Compiling problem file ..."
	$(COMPILER) $(PARAM) $(MAIN) $(USR) $(PHYSICS) $(HEADERS) $(FLAGS) -o $(EXEC)  
	@echo ""
	@echo "aztekas compiled successfully"
	
clean:
	rm -f $(EXEC)

###############################################
=======
###############################################
# Include a global Makefile in the /Src/ 
# directory, which include some definitions
# of parameters, inclusions of other makefiles
# and compilation orders
###############################################
include $(AZTPATH)/Src/Makefile
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
