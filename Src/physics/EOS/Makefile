##############################################
####### Makefile for Equation of State #######
##############################################
#---------------------------------------------
# Files are automatically included in
# compilation when the EOS variable is 
# specified. But, you can add them here.
#---------------------------------------------

ifeq ($(EOS),User)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(EOS),Ideal)
        PARAM += -DEOS=IDEAL
endif
ifeq ($(EOS),Dust)
        PARAM += -DEOS=DUST
endif
ifeq ($(EOS),Stiff)
        PARAM += -DEOS=STIFF	
endif
ifeq ($(EOS),Ryu)
        PARAM += -DEOS=RYU
endif

ifneq ($(EOS),User)
   include $(AZTPATH)/Src/physics/EOS/$(EOS)/Makefile
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
