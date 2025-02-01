##################################################
##             AZTEKAS USER MAKEFILE             #
##################################################

##################################################
# General Configuration
##################################################

# Name of the executable file to generate
EXECUTABLE = aztekas.run

# Directory where object files are stored
OBJDIR = ./obj

# Collect all object files in the object directory
OBJS := $(wildcard $(OBJDIR)/*.o)

##################################################
# OpenMP Configuration
##################################################
# Set the number of threads for OpenMP.
# To find the number of physical cores on your
# system, run:
# $ grep '^core id' /proc/cpuinfo | sort -u | wc -l
OMP_NUM = 4

##################################################
# Physics Configuration
##################################################
# Physics options available:
# - HD: Hydro (Euler inviscid equations)
# - RHD: Relativistic Hydro
PHY = HD

# Enable or disable specific physics modules
HYDRO = True
TOV   = False

##################################################
# Dimensionality
##################################################
# Set the number of dimensions:
# - 1: One-dimensional
# - 2: Two-dimensional
# - 3: Three-dimensional
# - 4: Two and a half dimensions
DIM = 1

##################################################
# Metric Selection
##################################################
# Choose the metric for relativistic physics:
# - Minkowski, Schwarzschild, Eddington-Finkelstein, Boyer-Lindquist, Kerr-Schild
# This setting is applicable only if PHY == RHD or RMHD
METRIC = Minkowski

##################################################
# Equation of State
##################################################
# Available options:
# - User
# - Ideal
# - Dust
# - Stiff
EOS = Ideal

##################################################
# Coordinate System
##################################################
# Select the coordinate system:
# - Cartesian
# - Cylindrical
# - Spherical
COORD = Spherical

##################################################
# Source Terms
##################################################
# Enable or disable user-defined source terms
USER_SOURCE = False

##################################################
# Integration Method
##################################################
# Specify the integration method:
# - Standard
INT = Standard

##################################################
# Local Source Files
##################################################
# Define local source files and their corresponding object files
LOCAL_SRCS = initial.cpp \
		user_boundaries.cpp \
		user_input.cpp \
		user_sources.cpp
LOCAL_OBJS = $(LOCAL_SRCS:%.cpp=$(OBJDIR)/%.o)

# Rule to compile .cpp files into .o files
$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

##################################################
# Global Include
##################################################
# Include a global Makefile from the /src/ directory
# This file contains additional parameters, definitions, and rules.
# Only modify this if you are familiar with the project's structure.
include $(AZTEKAS_PATH)/Makefile

##################################################
# Utility Targets
##################################################

.PHONY: init clean pre-commit commit docs help

# Install dependencies and activate the virtual environment with Poetry
init:

# Clean target: Add commands to remove compiled objects and temporary files
clean:

# Commit target: Placeholder for custom commit actions
commit:

# Build the last executable
last: last_$(EXECUTABLE)

last_$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@

# Pre-commit checks
pre-commit:

## Display available commands and descriptions
help:
