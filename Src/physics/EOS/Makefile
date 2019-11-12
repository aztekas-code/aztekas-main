##############################################
####### Makefile for Equation of State #######
##############################################
#---------------------------------------------
# Files are automatically included in
# compilation when the EOS variable is 
# specified. But, you can add them here.
#---------------------------------------------

ifneq ($(EOS),User)
   PHYSICS += $(AZTPATH)/Src/physics/EOS/$(EOS)/eos.c
endif

ifeq ($(EOS),User)
   $(info #######################)
   $(info # USER-DEFINED-EoS#)
   $(info #######################)
   $(info         )
   $(info In order to use a User Defined Equation of State (EoS) in aztekas,)
   $(info you will have to create the file user_eos.c)
   $(info inside your current directory.)
   $(info       )

   USR += user_eos.c

   ifeq (,$(wildcard ./user_eos.c))
      $(info No user_eos.c file)
      exit:
	      @echo "Exit. Compilation failed."
	      @echo " "
   endif
endif
