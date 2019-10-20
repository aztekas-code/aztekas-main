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

AZTPATH = $(AZTEKAS_PATH)

###############################################
# Physics (hd, mhd, rhd, rmhd)
###############################################

PHY = rhd

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

DIM = 2

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = cylindrical

###############################################
# User Source Terms (true,false)
###############################################

USER_SOURCE = false

###############################################
# Integration method
###############################################

INT = standard

###############################################

USR = boundaries.c \
		initial.c \
		user_input.c 

include $(AZTPATH)/Src/Makefile

###############################################
# Compilation
###############################################

FLAGS = -Ofast -lm -fopenmp
PARAM += -DOMP_NUM=6

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
