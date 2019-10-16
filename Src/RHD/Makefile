PHYSICS = $(AZTPATH)/Src/RHD/fvector.c \
			 $(AZTPATH)/Src/RHD/gvector.c \
		 	 $(AZTPATH)/Src/RHD/hvector.c \
			 $(AZTPATH)/Src/RHD/svector.c \
			 $(AZTPATH)/Src/RHD/qvector.c \
			 $(AZTPATH)/Src/RHD/q2uvector.c \
			 $(AZTPATH)/Src/RHD/u2qvector.c \
			 $(AZTPATH)/Src/RHD/surface.c \
		    $(AZTPATH)/Src/RHD/$(METRIC)/metric.c

ifeq ($(INT),pvrs)
   PHYSICS += $(AZTPATH)/Src/RHD/amatrix.c
endif
