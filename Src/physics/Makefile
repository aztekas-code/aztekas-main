###############################################
####### Makefile for Physics in AZTEKAS #######
###############################################
#----------------------------------------------
# All the physics you want to add is included
# here. But first, you have to specified in 
# the local Makefile, defining the 
# corresponding variable that will ensure the
# inclusion or certain kind of physics.
#----------------------------------------------

###############################################
# HYDRODYNAMICS 
###############################################

PARAM += -DPHYSICS=$(PHY)
include $(AZTPATH)/Src/physics/$(PHY)/Makefile

###############################################
# Equation of state
###############################################

ifeq ($(EOS),ideal)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(EOS),dust)
        PARAM += -DEOS=DUST
endif
ifeq ($(EOS),stiff)
        PARAM += -DEOS=STIFF	
endif

include $(AZTPATH)/Src/physics/EOS/Makefile
