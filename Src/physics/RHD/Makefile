ifeq ($(METRIC),User)
   PARAM += -DMETRIC=USER
endif
ifeq ($(METRIC),Minkowski)
   PARAM += -DMETRIC=MINK
endif
ifeq ($(METRIC),Schwarzschild)
   PARAM += -DMETRIC=SCHW
endif
ifeq ($(METRIC),Eddington-Finkelstein)
   PARAM += -DMETRIC=EF
endif
ifeq ($(METRIC),Boyer-Lindquist)
   PARAM += -DMETRIC=BL
endif
ifeq ($(METRIC),Kerr-Schild)
   PARAM += -DMETRIC=KS
endif

PHYSICS = $(AZTPATH)/Src/physics/RHD/fvector.c \
			 $(AZTPATH)/Src/physics/RHD/gvector.c \
		 	 $(AZTPATH)/Src/physics/RHD/hvector.c \
			 $(AZTPATH)/Src/physics/RHD/svector.c \
			 $(AZTPATH)/Src/physics/RHD/qvector.c \
			 $(AZTPATH)/Src/physics/RHD/q2uvector.c \
			 $(AZTPATH)/Src/physics/RHD/u2qvector.c \
			 $(AZTPATH)/Src/physics/RHD/surface.c \
		    $(AZTPATH)/Src/physics/RHD/$(METRIC)/metric.c

ifeq ($(INT),pvrs)
   PHYSICS += $(AZTPATH)/Src/RHD/amatrix.c
endif
