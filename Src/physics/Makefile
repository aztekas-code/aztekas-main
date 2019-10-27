#####################################

PARAM += -DPHYSICS=$(PHY)
include $(AZTPATH)/Src/physics/$(PHY)/Makefile

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

include $(AZTPATH)/Src/physics/EOS/Makefile
