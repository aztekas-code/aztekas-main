MAIN += $(AZTPATH)/Src/common/boundaries.c \
		  $(AZTPATH)/Src/common/default_boundaries.c \
        $(AZTPATH)/Src/common/freq_output.c \
        $(AZTPATH)/Src/common/initialize.c \
        $(AZTPATH)/Src/common/input.c \
        $(AZTPATH)/Src/common/messages.c \
        $(AZTPATH)/Src/common/output.c \
        $(AZTPATH)/Src/common/restart.c \
        $(AZTPATH)/Src/common/termination.c \
        $(AZTPATH)/Src/common/tools.c

ifeq ($(HYDRO),True)
   MAIN += $(AZTPATH)/Src/common/print_time_values.c
endif
ifeq ($(TOV),True)
   MAIN += $(AZTPATH)/Src/common/print_values.c
endif
