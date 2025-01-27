##################################################
#####          MAKEFILE FOR AZTEKAS          #####
##################################################
#-------------------------------------------------
# From the list below, please activate/deactivate 
# the options that apply to your run. If you 
# modify any of these options, make sure that you 
# recompile the whole code by typing 
# $ make clean && make 
#-------------------------------------------------

##################################################
# PATH TO AZTEKAS
##################################################

EXECUTABLE = aztekas.run
OBJDIR = ./obj
OBJS := $(wildcard $(OBJDIR)/*.o)

##################################################
# SET NUMBER OF THREADS USED BY OMP
##################################################
# You can find the number of physical cores 
# using
# $ grep '^core id' /proc/cpuinfo |sort -u|wc -l
##################################################

OMP_NUM = 4

##################################################
# Physics (HD)
#    - HYDRO (Euler inviscid equations)
#    - Lane-Endem?
# Physics (RHD)
#    - HYDRO (Euler inviscid equations)
#    - TOV (Tolman-Oppenheimer-Volkoff, RHD only)
##################################################

PHY = HD

HYDRO = True
TOV   = False

##################################################
# Dimension (1, 2 or 3. Enter 4 for 2.5)
##################################################

DIM = 1

##################################################
# Metric (User, Minkowski, Schwarzschild, 
# Eddington-Finkelstein, Boyer-Lindquist or
# Kerr-Schild)
# ONLY IF PHY == RHD or RMHD
##################################################

METRIC = Minkowski

##################################################
# Equation of State (User, Ideal, Dust, Stiff)
##################################################

EOS = Ideal

##################################################
# Coordinates (Cartesian, Cylindrical, Spherical)
##################################################

COORD = Spherical

##################################################
# User Source Terms (True,False)
##################################################

USER_SOURCE = False

##################################################
# Integration method
##################################################

INT = Standard

##################################################
# All user *.cpp are written here
##################################################

LOCAL_SRCS = initial.cpp
LOCAL_OBJS = $(LOCAL_SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

##################################################
# Include a global Makefile in the /src/ 
# directory, which has some definitions
# of parameters, inclusions of other makefiles
# and compilation orders. Change this Makefile,
# and recursive Makefiles only if you know what
# your doing.
##################################################
include $(AZTEKAS_PATH)/Makefile

clean:

commit:

pre-commit:

previous: previous_$(EXECUTABLE)

previous_$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
