#####################################

ifeq ($(PHY),hd)
      PARAM += -DPHYSICS=HD
		include $(AZTPATH)/Src/HD/Makefile
endif
ifeq ($(PHY),rhd)
      PARAM += -DPHYSICS=RHD
		include $(AZTPATH)/Src/RHD/Makefile
endif

#####################################

ifeq ($(EOS),ideal)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(EOS),dust)
        PARAM += -DEOS=DUST
endif
ifeq ($(EOS),stiff)
        PARAM += -DEOS=STIFF	
endif

include $(AZTPATH)/Src/EOS/Makefile

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

ifeq ($(USER_SOURCE),true)
   PARAM += -DUSER_SOURCE_TERMS=TRUE
   USR += user_sources.c
else
   PARAM += -DUSER_SOURCE_TERMS=FALSE
endif

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

PARAM += -Dgc=3

#####################################

MAIN = $(AZTPATH)/Src/alloc.c \
       $(AZTPATH)/Src/auxfunc.c \
       $(AZTPATH)/Src/bound_cond.c \
       $(AZTPATH)/Src/flux.c \
       $(AZTPATH)/Src/input.c \
       $(AZTPATH)/Src/integration$(INT).c \
       $(AZTPATH)/Src/limiters.c \
		 $(AZTPATH)/Src/main.c \
		 $(AZTPATH)/Src/mesh.c \
		 $(AZTPATH)/Src/output.c \
		 $(AZTPATH)/Src/restart.c \
		 $(AZTPATH)/Src/runge-kutta.c \
		 $(AZTPATH)/Src/timestep.c \
		 $(AZTPATH)/Src/vectors.c 
