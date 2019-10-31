#####################################
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
	PARAM += -Deq=3
	PARAM += -Dgraf=1
endif
ifeq ($(DIM),2)
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif
ifeq ($(DIM),4)
	PARAM += -Deq=5
	PARAM += -Dgraf=2
endif

PARAM += -DDIM=$(DIM)

#####################################
# COORDINATES
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
# USER SOURCE TERMS 
#####################################

ifeq ($(USER_SOURCE),true)
   $(info ########################)
   $(info # USER-DEFINED-SOURCES #)
   $(info ########################)
   $(info In order to use a User Defined Sources in aztekas,)
   $(info you will have to create the file user_sources.c)
   $(info inside your current directory. Remember that the)
   $(info user_sources.c file must have one void function)
   $(info void User_Source_Terms(double *s, double *u, gauge_ *local_grid))

   PARAM += -DUSER_SOURCE_TERMS=TRUE
   USR += user_sources.c

   ifeq (,$(wildcard ./user_sources.c))
      $(info       )
      $(info No user_sources.c file)
      exit:
	   	@echo "Exit failed compilation"
endif

else
   PARAM += -DUSER_SOURCE_TERMS=FALSE
endif

#####################################
# INTEGRATION 
#####################################

ifeq ($(INT),standard)
   PARAM += -DINTEGRATION=STANDARD
   INT = 
endif
ifeq ($(INT),pvrs)
   PARAM += -DINTEGRATION=PVRS
   INT = -pvrs
endif

#####################################
# GHOST CELLS
#####################################

PARAM += -Dgc=3

#####################################

#####################################
# Include recursive Makefiles
#####################################

MAIN = $(AZTPATH)/Src/main.c

include $(AZTPATH)/Src/common/Makefile
include $(AZTPATH)/Src/fluxes/Makefile
include $(AZTPATH)/Src/grid/Makefile
include $(AZTPATH)/Src/include/Makefile
include $(AZTPATH)/Src/integration/Makefile
include $(AZTPATH)/Src/physics/Makefile

###############################################
# Compilation
###############################################

THREADS = $(shell grep '^core id' /proc/cpuinfo |sort -u|wc -l)

FLAGS = -Ofast -lm
ifneq ($(OMP_NUM),1)
   FLAGS += -fopenmp
endif

COMPILER = gcc 
AZT_HEAD = $(AZTPATH)/Src/include
HEADERS = -I$(AZT_HEAD) -I.
EXEC = aztekas

aztekas: 
	@echo ""
	@echo "###############################################"
	@echo "##### AZTEKAS code Makefile compilation #######"
	@echo "###############################################"
	@echo ""
	@echo "Compiling files..."
	@echo ""
	$(COMPILER) $(PARAM) $(MAIN) $(USR) $(PHYSICS) $(HEADERS) $(FLAGS) -o $(EXEC)  
	@echo ""
	@echo "# AZTEKAS will run with" $(OMP_NUM)" of "$(THREADS) "cores available to use with OMP."
	@echo "# You can modify the OMP_NUM variable inside the Makefile."
	@echo ""
	@echo "###############################################"
	@echo "##### AZTEKAS code compiled successfully ######"
	@echo "###############################################"
	@echo ""
	
clean:
	rm -f $(EXEC)

###############################################
