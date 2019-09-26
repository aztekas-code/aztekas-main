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

AZTPATH = ../../../

###############################################
# Physics (hd, mhd, rhd, rmhd)
###############################################

PHY = hd

########################################################################################
# Metric (Minkowski, Schwarzschild, Eddington-Finkelstein, Boyer-Lindquist, Kerr-Schild)
########################################################################################

## ONLY IF PHY == rhd or rmhd

METRIC = Minkowski

###############################################
# Equation of State (ideal, dust, stiff)
###############################################

eos = ideal

###############################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
###############################################

DIM = 1

###############################################
# Coordinates (cartesian,cylindrical,spherical)
###############################################

COORD = cartesian

###############################################
# Integration method
###############################################

INT = standard

###############################################
# Ghost cells
###############################################

PARAM += -Dgc=3

###############################################

ifeq ($(PHY),hd)
        PARAM += -DPHYSICS=HD
endif
ifeq ($(PHY),rhd)
        PARAM += -DPHYSICS=RHD
endif

#####################################

ifeq ($(eos),ideal)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(eos),dust)
        PARAM += -DEOS=DUST
endif
ifeq ($(eos),stiff)
        PARAM += -DEOS=STIFF	
endif

#####################################

ifeq ($(DIM),1)
	PARAM += -DDIM=1
	PARAM += -Deq=3
	PARAM += -Dgraf=1
endif
ifeq ($(DIM),2)
	PARAM += -DDIM=2
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
	PARAM += -DDIM=3
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif
ifeq ($(DIM),4)
	PARAM += -DDIM=4
	PARAM += -Deq=5
	PARAM += -Dgraf=2
endif

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

ifeq ($(INT),standard)
	PARAM += -Dintegration=0
	INT = 
endif
ifeq ($(INT),pvrs)
	PARAM += -Dintegration=1
	INT = -pvrs
endif

#####################################

#####################################
# Compilation
#####################################

ifeq ($(PHY),hd)
SOURCES = $(AZTPATH)/Src/main.c \
	  $(AZTPATH)/Src/alloc.c \
	  $(AZTPATH)/Src/mesh.c \
	  ./initial.c \
	  $(AZTPATH)/Src/timestep.c \
	  $(AZTPATH)/Src/input.c \
	  $(AZTPATH)/Src/output.c \
	  $(AZTPATH)/Src/integration$(INT).c \
	  $(AZTPATH)/Src/flux.c \
	  $(AZTPATH)/Src/runge-kutta.c \
	  $(AZTPATH)/Src/limiters.c \
	  $(AZTPATH)/Src/vectors.c \
	  $(AZTPATH)/Src/restart.c \
	  $(AZTPATH)/Src/auxfunc.c \
	  ./boundaries.c \
	  $(AZTPATH)/Src/bound_cond.c \
	  $(AZTPATH)/Src/HD/q2uvector.c \
	  $(AZTPATH)/Src/HD/u2qvector.c \
	  $(AZTPATH)/Src/HD/qvector.c \
	  $(AZTPATH)/Src/HD/fvector.c \
	  $(AZTPATH)/Src/HD/gvector.c \
	  $(AZTPATH)/Src/HD/hvector.c \
	  $(AZTPATH)/Src/HD/svector.c \
	  $(AZTPATH)/Src/HD/surface.c \
	  $(AZTPATH)/Src/EOS/eos.c \
	  ./user_sources.c \
	  ./user_input.c
endif 

ifeq ($(PHY),rhd)
SOURCES = $(AZTPATH)/Src/main.c \
	  $(AZTPATH)/Src/alloc.c \
	  $(AZTPATH)/Src/mesh.c \
	  ./initial.c \
	  $(AZTPATH)/Src/timestep.c \
	  $(AZTPATH)/Src/input.c \
	  $(AZTPATH)/Src/output.c \
	  $(AZTPATH)/Src/integration$(INT).c \
	  $(AZTPATH)/Src/flux.c \
	  $(AZTPATH)/Src/runge-kutta.c \
	  $(AZTPATH)/Src/limiters.c \
	  $(AZTPATH)/Src/vectors.c \
	  $(AZTPATH)/Src/restart.c \
	  $(AZTPATH)/Src/auxfunc.c \
	  ./boundaries.c \
	  $(AZTPATH)/Src/bound_cond.c \
	  $(AZTPATH)/Src/RHD/q2uvector.c \
	  $(AZTPATH)/Src/RHD/u2qvector.c \
	  $(AZTPATH)/Src/RHD/qvector.c \
	  $(AZTPATH)/Src/RHD/fvector.c \
	  $(AZTPATH)/Src/RHD/gvector.c \
	  $(AZTPATH)/Src/RHD/hvector.c \
	  $(AZTPATH)/Src/RHD/svector.c \
	  $(AZTPATH)/Src/RHD/surface.c \
	  $(AZTPATH)/Src/RHD/$(METRIC)/metric.c \
	  $(AZTPATH)/Src/EOS/eos.c \
	  ./user_sources.c \
	  ./user_input.c
endif

#OBJS = $(SOURCES:.c=.o)

FLAGS = -Ofast -lm

COMPILER = gcc
AZT_HEAD = $(AZTPATH)/Src/Headers
CFLAGS = -I$(AZT_HEAD) -I.
EXEC = aztekas

#####################################

$(EXEC): $(SOURCES)
	@echo ""
	@echo "Compiling problem file ..."
	$(COMPILER) $(PARAM) -fopenmp $(SOURCES) $(FLAGS) $(CFLAGS) -o $(EXEC)  
	@echo "$(bla)"
	@echo "aztekas compiled successfully"
	
clean:
	rm -f $(EXEC)
	
