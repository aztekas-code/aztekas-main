#####################################
# GLOBAL MAKEFILE OF AZTEKAS
#####################################

#####################################
# OMP PARALELIZATION
#####################################

PARAM += -DOMP_NUM=$(OMP_NUM)

#####################################
# PHYSICS
#####################################

ifeq ($(PHY),HD)
   METRIC=-
endif

#####################################
# IMPLEMENTED SYSTEMS OF EQUATIONS 
#####################################

ifeq ($(HYDRO),True)
   PARAM += -DHYDRO=TRUE
   PARAM += -DHYPERBOLIC=TRUE
else
   HYDRO=False
   PARAM += -DHYDRO=FALSE
   PARAM += -DHYPERBOLIC=FALSE
endif
ifeq ($(TOV),True)
   PARAM += -DTOV=TRUE
   PARAM += -DODE=TRUE
else
   TOV=False
   PARAM += -DTOV=FALSE
   PARAM += -DODE=FALSE
endif

ifeq ($(HYDRO),False)
   ifeq ($(TOV),False)
      $(info #####################)
      $(info # COMPILATION ERROR #)
      $(info #####################)
      $(info )
      $(info Not system of equations to solve specified)
      $(info Include in your local Makefile at least one of the possible implemented problems:)
      $(info - Euler inviscid equations    -> HYDRO=True)
      $(info - Tolman-Oppenheimer-Volkoff  -> TOV=True)
      $(info )
      $(error Failed compilation)
   endif
endif

#####################################
# DIMENSIONS
#####################################

ifeq ($(DIM),1)
   PARAM += -Deq=3
   DIMENSION=1
endif
ifeq ($(DIM),2)
   PARAM += -Deq=4
   DIMENSION=2
endif
ifeq ($(DIM),3)
   PARAM += -Deq=5
   DIMENSION=3
endif
ifeq ($(DIM),4)
   PARAM += -Deq=5
   DIMENSION="2.5"
endif

PARAM += -DDIM=$(DIM)

#####################################
# COORDINATES
#####################################

ifeq ($(COORD),Cartesian)
	PARAM += -DCOORDINATES=CARTESIAN
endif
ifeq ($(COORD),Cylindrical)
	PARAM += -DCOORDINATES=CYLINDRICAL
endif
ifeq ($(COORD),Spherical)
	PARAM += -DCOORDINATES=SPHERICAL
endif

#####################################
# USER SOURCE TERMS 
#####################################

ifeq ($(USER_SOURCE),True)
   $(info ########################)
   $(info # USER-DEFINED-SOURCES #)
   $(info ########################)
   $(info In order to use a User Defined Sources in aztekas,)
   $(info you will have to create the file user_sources.c)
   $(info inside your current directory. Remember that the)
   $(info user_sources.c file must have one void function)
   $(info void User_Source_Terms(double *s, double *u, gauge_ *local_grid))
   $(info       )

   PARAM += -DUSER_SOURCE_TERMS=TRUE
   USR += user_sources.c

   ifeq (,$(wildcard ./user_sources.c))
      $(info No user_sources.c file)
      exit:
	   	@echo "Exit. Compilation failed."
			@echo " "

	$(info user_sources.c file does exists!)
   $(info       )
endif

else
   PARAM += -DUSER_SOURCE_TERMS=FALSE
endif

#####################################
# INTEGRATION 
#####################################

ifneq ($(INT),$(filter $(INT),Standard PVRS))
   exit:
	   @echo "Wrong integration parameter."
		@echo "Possible options:"
		@echo "   - Standard"
		@echo "   - PVRS"
	   @echo "Exit. Compilation failed."
	   @echo " "
endif

ifeq ($(INT),Standard)
   PARAM += -DINTEGRATION=STANDARD
   INT = 
   INTEGRATION="Standard FVM"
endif
ifeq ($(INT),PVRS)
   PARAM += -DINTEGRATION=PVRS
   INT = -pvrs
   INTEGRATION=PVRS
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
include $(AZTPATH)/Src/integration/Makefile
include $(AZTPATH)/Src/physics/Makefile

###############################################
# Compilation
###############################################

THREADS = $(shell grep '^core id' /proc/cpuinfo |sort -u|wc -l)

FLAGS = -lm -fcommon -std=gnu99
ifneq ($(OMP_NUM),1)
   FLAGS += -fopenmp
endif

COMPILER = gcc
AZT_HEAD = $(AZTPATH)/Src/include
HEADERS = -I$(AZT_HEAD) -I.
EXEC = aztekas


all: 
	@echo ""
	@echo "###########################################################"
	@echo "########### AZTEKAS code Makefile compilation #############"
	@echo "###########################################################"
	@echo ""
	@echo "Compiling files..."
	@echo ""
	@-$(COMPILER) $(PARAM) $(MAIN) $(USR) $(PHYSICS) $(HEADERS) $(FLAGS) -o $(EXEC) && echo "\033[0;32m=== Compilation for $(EXEC) SUCCEED ===\033[0m" || echo "\033[1;31m=== Compilation $(EXEC) FAILED ===\033[0m";
	@echo ""
	@echo "###########################################################"
	@echo "############ AZTEKAS code Makefile parameters #############"
	@echo "###########################################################"
	@echo "" 
	@echo " AZTEKAS code path = "$(AZTPATH)
	@echo " Number of OMP cores = "$(OMP_NUM)
	@echo " Physics = "$(PHY)
	@echo " Equation of State (EoS) = "$(EOS)
	@echo " Systems of equations"
	@echo "    - Euler inviscid = "$(HYDRO)
	@echo "    - Tolman-Oppenheimer-Volkoff = "$(TOV)
	@echo " Integration method = "$(INTEGRATION)
	@echo " Metric = "$(METRIC)
	@echo " Dimensions = "$(DIMENSION)
	@echo " Coordinate system = "$(COORD)
	@echo " User sources = " $(USER_SOURCE)
	@echo " User files = "$(USR) 
	@echo ""
	@echo "###########################################################"
	@echo "" 
	
clean:
	rm -f $(EXEC)

###############################################
