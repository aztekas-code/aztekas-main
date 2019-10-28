#####################################
<<<<<<< HEAD

ifeq ($(DIM),1)
	PARAM += -DDIM=1
=======
# GLOBAL MAKEFILE OF AZTEKAS
#####################################

#####################################
# OMP PARALELIZATION
#####################################

PARAM += -DOMP_NUM=$(OMP_NUM)

#####################################
# DIMENSIONS
#####################################

ifeq ($(DIM),1)
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
	PARAM += -Deq=3
	PARAM += -Dgraf=1
endif
ifeq ($(DIM),2)
<<<<<<< HEAD
	PARAM += -DDIM=2
=======
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
<<<<<<< HEAD
	PARAM += -DDIM=3
=======
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif
ifeq ($(DIM),4)
<<<<<<< HEAD
	PARAM += -DDIM=4
=======
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
	PARAM += -Deq=5
	PARAM += -Dgraf=2
endif

<<<<<<< HEAD
=======
PARAM += -DDIM=$(DIM)

#####################################
# COORDINATES
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
#####################################

ifeq ($(COORD),cartesian)
	PARAM += -DCOORDINATES=CARTESIAN
endif
ifeq ($(COORD),cylindrical)
	PARAM += -DCOORDINATES=CYLINDRICAL
endif
ifeq ($(COORD),spherical)
	PARAM += -DCOORDINATES=SPHERICAL
endif

#####################################
<<<<<<< HEAD
=======
# USER SOURCE TERMS 
#####################################
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

ifeq ($(USER_SOURCE),true)
   PARAM += -DUSER_SOURCE_TERMS=TRUE
   USR += user_sources.c
else
   PARAM += -DUSER_SOURCE_TERMS=FALSE
endif

#####################################
<<<<<<< HEAD
=======
# INTEGRATION 
#####################################
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

ifeq ($(INT),standard)
	PARAM += -DINTEGRATION=STANDARD
	INT = 
endif
ifeq ($(INT),pvrs)
	PARAM += -DINTEGRATION=PVRS
	INT = -pvrs
endif

#####################################

PARAM += -Dgc=3

#####################################

<<<<<<< HEAD
MAIN = $(AZTPATH)/Src/main.c \
=======
MAIN = $(AZTPATH)/Src/main.c

>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907

include $(AZTPATH)/Src/common/Makefile
include $(AZTPATH)/Src/fluxes/Makefile
include $(AZTPATH)/Src/grid/Makefile
include $(AZTPATH)/Src/integration/Makefile
include $(AZTPATH)/Src/physics/Makefile
<<<<<<< HEAD
=======

###############################################
# Compilation
###############################################

FLAGS = -Ofast -lm 
ifneq ($(OMP_NUM),1)
   FLAGS += -fopenmp
endif

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
>>>>>>> 12b3acd607466560c2bebb7b61677f23252c7907
