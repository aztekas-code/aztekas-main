MAIN += $(AZTPATH)/Src/integration/equation-solver.c \
	     $(AZTPATH)/Src/integration/runge-kutta.c

ifeq ($(HYDRO),True)
   include  $(AZTPATH)/Src/integration/HYPERBOLIC/Makefile
endif
ifeq ($(TOV),True)
   include  $(AZTPATH)/Src/integration/ODE/Makefile
endif
