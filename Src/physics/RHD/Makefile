
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

PHYSICS = $(AZTPATH)/Src/physics/RHD/cons_q.c \
			 $(AZTPATH)/Src/physics/RHD/cons2prim.c \
			 $(AZTPATH)/Src/physics/RHD/flux_f.c \
			 $(AZTPATH)/Src/physics/RHD/flux_g.c \
		 	 $(AZTPATH)/Src/physics/RHD/flux_h.c \
			 $(AZTPATH)/Src/physics/RHD/prim2cons.c \
			 $(AZTPATH)/Src/physics/RHD/sources.c \
			 $(AZTPATH)/Src/physics/RHD/surface_volume.c

ifneq ($(METRIC),User)
   PHYSICS += $(AZTPATH)/Src/physics/RHD/$(METRIC)/metric.c
endif

ifeq ($(METRIC),User)
   $(info #######################)
   $(info # USER-DEFINED-METRIC #)
   $(info #######################)
   $(info         )
   $(info In order to use a User Defined Metric in aztekas,)
   $(info you will have to create the file user_metric.c)
   $(info inside your current directory. Remember that the)
   $(info user_metric.c file must have two void functions)
   $(info void Get_Metric_Components(gauge_ *local_grid))
   $(info and)
   $(info void Gauge_Derivatives(der_gauge_ *der, gauge_ *local_grid))
   $(info in which you define the components of the metric in terms)
   $(info of the 3+1 formalism [lapse, shift and 3D-metric])
   $(info and their derivatives, respectively.)
   $(info       )

   PHYSICS += user_metric.c

   ifeq (,$(wildcard ./user_metric.c))
      $(info No user_metric.c file)
      exit:
	      @echo "Exit. Compilation failed."
	      @echo " "
   endif
endif




ifeq ($(INT),pvrs)
   PHYSICS += $(AZTPATH)/Src/RHD/amatrix.c
endif
