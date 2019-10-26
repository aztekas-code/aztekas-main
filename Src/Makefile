#####################################

PARAM += -DOMP_NUM=$(OMP_NUM)

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

MAIN = $(AZTPATH)/Src/main.c \

include $(AZTPATH)/Src/common/Makefile
include $(AZTPATH)/Src/fluxes/Makefile
include $(AZTPATH)/Src/grid/Makefile
include $(AZTPATH)/Src/integration/Makefile
include $(AZTPATH)/Src/physics/Makefile
