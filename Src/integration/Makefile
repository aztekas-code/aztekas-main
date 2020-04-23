ifeq ($(PHY),HD) 
   MAIN += $(AZTPATH)/Src/integration/equation-solver.c \
	   	  $(AZTPATH)/Src/integration/integration$(INT).c \
		     $(AZTPATH)/Src/integration/method-of-lines.c \
		     $(AZTPATH)/Src/integration/runge-kutta.c \
		     $(AZTPATH)/Src/integration/timestep.c
endif
ifeq ($(PHY),RHD) 
   MAIN += $(AZTPATH)/Src/integration/equation-solver.c \
	   	  $(AZTPATH)/Src/integration/integration$(INT).c \
		     $(AZTPATH)/Src/integration/method-of-lines.c \
		     $(AZTPATH)/Src/integration/runge-kutta.c \
		     $(AZTPATH)/Src/integration/timestep.c
endif
ifeq ($(PHY),TOV)
   MAIN += $(AZTPATH)/Src/integration/equaton-solver.c \
			  $(AZTPATH)/Src/integration/integration.c \
		     $(AZTPATH)/Src/integration/runge-kutta.c
endif
