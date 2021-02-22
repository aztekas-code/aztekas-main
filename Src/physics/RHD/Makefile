######################################################
###### Makefile for Relativistic Hydrodynamics #######
######################################################

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

   USR += user_metric.c

   ifeq (,$(wildcard ./user_metric.c))
      $(info No user_metric.c file)
      exit:
	      @echo "Exit. Compilation failed."
	      @echo " "
   endif
endif

# Systems of equations

ifeq ($(HYDRO),True)
	include $(AZTPATH)/Src/physics/RHD/HYDRO/Makefile
endif

ifeq ($(TOV),True)
	include $(AZTPATH)/Src/physics/RHD/TOV/Makefile
endif
