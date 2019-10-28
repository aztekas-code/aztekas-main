#####################################
#### Define dimension parameters ####
#####################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------

###############################################
# PATH TO AZTEKAS
###############################################
<<<<<<< HEAD

AZTPATH = $(AZTEKAS_PATH)

=======

AZTPATH = $(AZTEKAS_PATH)

###############################################
# SET NUMBER OF THREADS USED BY OMP
###############################################

OMP_NUM = 6

>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
###############################################
# Physics (hd, mhd, rhd, rmhd)
###############################################

PHY = hd

########################################################################################
# Metric (Minkowski, Schwarzschild, Eddington-Finkelstein, Boyer-Lindquist, Kerr-Schild)
########################################################################################

## ONLY IF PHY == rhd or rmhd

METRIC = Minkowski
<<<<<<< HEAD

###############################################
# Equation of State (ideal, dust, stiff)
###############################################

EOS = ideal

###############################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 2

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = spherical

###############################################
# User Source Terms (true,false)
###############################################

=======

###############################################
# Equation of State (ideal, dust, stiff)
###############################################

EOS = ideal

###############################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 2

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = spherical

###############################################
# User Source Terms (true,false)
###############################################

>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
USER_SOURCE = true

###############################################
# Integration method
###############################################

INT = standard

###############################################
<<<<<<< HEAD

USR = boundaries.c \
		initial.c \
		user_input.c 

include $(AZTPATH)/Src/Makefile

###############################################
# Compilation
###############################################

FLAGS = -Ofast -lm -fopenmp
PARAM += -DOMP_NUM=2
=======

USR = boundaries.c \
		initial.c \
		user_input.c 

include $(AZTPATH)/Src/Makefile

###############################################
# Compilation
###############################################

FLAGS = -Ofast -lm
ifneq ($(OMP_NUM),1)
   FLAGS += -fopenmp
endif
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

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
