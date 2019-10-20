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
