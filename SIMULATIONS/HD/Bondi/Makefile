######################################
#### Define simulation parameters ####
######################################
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

PHY = hd
ifeq ($(PHY),hd)
        PARAM += -DPHYSICS=1
endif
ifeq ($(PHY),rhd)
        PARAM += -DPHYSICS=2
endif

#####################################

#####################################
# Dimension
#####################################

DIM = 1

ifeq ($(DIM),1)
	PARAM += -Ddim=1
	PARAM += -Deq=3
	PARAM += -Dgraf=1
endif
ifeq ($(DIM),2)
	PARAM += -Ddim=2
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),4)
	PARAM += -Ddim=2
	PARAM += -Deq=4
	PARAM += -Dgraf=2
endif
ifeq ($(DIM),3)
	PARAM += -Ddim=3
	PARAM += -Deq=5
	PARAM += -Dgraf=3
endif

#####################################

#####################################
# Coordinates (cart,cyl,sph)
#####################################

COORD = sph

ifeq ($(COORD),cart)
	PARAM += -Dalfa=0
endif
ifeq ($(COORD),cyl)
	PARAM += -Dalfa=1
endif
ifeq ($(COORD),sph)
	PARAM += -Dalfa=2
endif

#####################################

#####################################
# PATH TO AZTEKAS
#####################################
AZTPATH = /home/aaguayoo/DriveCiencias/Articles/AZTEKAS/aguayo-mendoza_aztekas/code
#####################################

#####################################
# Integration method
#####################################

#####################################

#####################################
# Compilation
#####################################

ifeq ($(PHY),hd)
SOURCES = $(AZTPATH)/Src/main.c \
	  $(AZTPATH)/Src/alloc.c \
	  $(AZTPATH)/Src/array.c \
	  $(AZTPATH)/Src/mesh.c \
	  ./initial.c \
	  $(AZTPATH)/Src/timestep.c \
	  $(AZTPATH)/Src/input.c \
	  $(AZTPATH)/Src/output.c \
	  $(AZTPATH)/Src/integration.c \
	  $(AZTPATH)/Src/flux.c \
	  $(AZTPATH)/Src/rk-clas-$(COORD).c \
	  $(AZTPATH)/Src/limiters.c \
	  $(AZTPATH)/Src/vectors.c \
	  $(AZTPATH)/Src/restart.c\
	  $(AZTPATH)/Src/auxfunc.c\
	  ./boundaries.c \
	  $(AZTPATH)/Src/bound_cond.c \
	  $(AZTPATH)/Src/HD-$(COORD)/amatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/dmmatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/dnmatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/domatrix.c \
	  $(AZTPATH)/Src/HD-$(COORD)/q2uvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/u2qvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/qvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/fvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/gvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)/hvector.c \
	  $(AZTPATH)/Src/HD-$(COORD)$(grav)/svector.c \
	  $(AZTPATH)/Src/HD-$(COORD)$(grav)/gauge.c \
	  ./extforce.c
endif 

ifeq ($(PHY),rhd)
SOURCES = $(AZTPATH)/Src/main.c \
	  $(AZTPATH)/Src/alloc.c \
	  $(AZTPATH)/Src/array.c \
	  $(AZTPATH)/Src/mesh.c \
	  ./initial.c \
	  $(AZTPATH)/Src/timestep.c \
	  $(AZTPATH)/Src/input.c \
	  $(AZTPATH)/Src/output.c \
	  $(AZTPATH)/Src/integration.c \
	  $(AZTPATH)/Src/flux.c \
	  $(AZTPATH)/Src/rk-rel.c \
	  $(AZTPATH)/Src/limiters.c \
	  $(AZTPATH)/Src/vectors.c \
	  $(AZTPATH)/Src/restart.c\
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
	  $(AZTPATH)/Src/RHD-$(COORD)/gauge.c  
	  ./extforce.c
endif

#OBJS = $(SOURCES:.c=.o)

FLAGS = -Ofast -lm

COMPILER = gcc
AZT_HEAD = /home/aaguayoo/DriveCiencias/Articles/AZTEKAS/aguayo-mendoza_aztekas/code/Src/Headers
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
