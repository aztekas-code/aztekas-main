###################################################
###### Makefile for Classical Hydrodynamics #######
###################################################

PHYSICS = $(AZTPATH)/Src/physics/HD/cons_q.c \
			 $(AZTPATH)/Src/physics/HD/cons2prim.c \
			 $(AZTPATH)/Src/physics/HD/flux_f.c \
			 $(AZTPATH)/Src/physics/HD/flux_g.c \
		 	 $(AZTPATH)/Src/physics/HD/flux_h.c \
			 $(AZTPATH)/Src/physics/HD/prim2cons.c \
			 $(AZTPATH)/Src/physics/HD/sources.c \
			 $(AZTPATH)/Src/physics/HD/surface_volume.c

ifeq ($(INT),pvrs)
   PHYSICS += $(AZTPATH)/Src/physics/HD/matrix_a.c
endif
