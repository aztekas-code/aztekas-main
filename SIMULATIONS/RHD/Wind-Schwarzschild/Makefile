#####################################
#### Define dimension parameters ####
#####################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------

#####################################
# Ghost cells
#####################################

PARAM += -Dgc=3

#####################################

#####################################
# Physics (rhd,hd,mhd,rmhd)
#####################################

PHY    = rhd

# Metric (Minkowski, Schwarzschild, Eddington-Finkelstein, Kerr-Schild)
METRIC = Schwarzschild

ifeq ($(PHY),hd)
        PARAM += -DPHYSICS=HD
endif
ifeq ($(PHY),rhd)
        PARAM += -DPHYSICS=RHD
endif

eos = ideal
ifeq ($(eos),ideal)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(eos),dust)
        PARAM += -DEOS=DUST
endif

#####################################

#####################################
# Dimension
#####################################

DIM = 2

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
ifeq ($(DIM),4)
	PARAM += -DDIM=4
	PARAM += -Deq=5
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
	PARAM += -DDIM=3
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif

#####################################

#####################################
# Coordinates (cart,cyl,sph)
#####################################

COORD = sph

ifeq ($(COORD),cart)
	PARAM += -DCOORDINATES=CARTESIAN
endif
ifeq ($(COORD),cyl)
	PARAM += -DCOORDINATES=CYLINDRICAL
endif
ifeq ($(COORD),sph)
	PARAM += -DCOORDINATES=SPHERICAL
endif

#####################################

#####################################
# PATH TO AZTEKAS
#####################################
AZTPATH = ../../../
#####################################

#####################################
# Integration method
#####################################

INT = standard

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
	  $(AZTPATH)/Src/HD/metric.c \
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
	
