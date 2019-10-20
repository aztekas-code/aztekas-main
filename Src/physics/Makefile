#####################################

ifeq ($(PHY),hd)
      PARAM += -DPHYSICS=HD
      include $(AZTPATH)/Src/physics/HD/Makefile
endif
ifeq ($(PHY),rhd)
      PARAM += -DPHYSICS=RHD
      include $(AZTPATH)/Src/physics/RHD/Makefile
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

include $(AZTPATH)/Src/physics/EOS/Makefile
