#####################################
#### define dimension parameters ####
#####################################
#----------------------------------------------------------------------
# From the list below, please activate/deactivate the options that     
# apply to your run. If you modify any of these options, make sure     
# that you recompile the whole code by typing "make clean; make" 
#----------------------------------------------------------------------

#####################################
# Tests
# 0. Custom problem
# 1. Riemann problem (1D, 2D, 3D)
# 2. Kelvin-Helmholtz (2D, 3D)
# 3. Jet (1D, 2D, 3D)
# 4. Spherical Accretion (1D, 2D, 3D)
# 5. Wind (2D, 3D)
#####################################

# Set test 
PARAM += -DTEST=1

#####################################
# Ghost cells
#####################################

PARAM += -Dgc=3

#####################################

#####################################
# Physics (rhd,hd,mhd,rmhd)
#####################################

PHY = hd
ifeq ($(PHY),hd)
        PARAM += -DPHYSICS=HD
endif
ifeq ($(PHY),rhd)
        PARAM += -DPHYSICS=RHD
endif

# -Grav if you want classical gravity
# nothing if you don't.
grav = #-grav

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

COORD = cart

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
	  $(AZTPATH)/Src/HD-$(COORD)/dmmatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/dnmatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/domatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/q2uvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/u2qvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/qvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/fvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/gvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/hvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/svector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/gauge.c \
	  $(AZTPATH)/Src/HD-$(COORD)/metric.c \
	  ./extforce.c \
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
	  $(AZTPATH)/Src/rk-rel.c \
	  $(AZTPATH)/Src/limiters.c \
	  $(AZTPATH)/Src/vectors.c \
	  $(AZTPATH)/Src/restart.c \
	  $(AZTPATH)/Src/auxfunc.c \
	  ./boundaries.c \
	  $(AZTPATH)/Src/bound_cond.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/amatrix.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/dmmatrix.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/dnmatrix.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/domatrix.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/q2uvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/u2qvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/qvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/fvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/gvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/hvector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/svector.c \
	  $(AZTPATH)/Src/RHD-$(COORD)/gauge.c \
	  ./extforce.c
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
	
