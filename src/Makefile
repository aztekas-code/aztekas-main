# src/Makefile
include $(AZTEKAS_PATH)/src/common/Makefile
include $(AZTEKAS_PATH)/src/grid/Makefile

# Si no se define desde afuera, será "."
OBJDIR ?= .

# Detectar el sistema operativo
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # Estamos en macOS
  CXX      = clang++
  CXXFLAGS = -Wall -std=c++17 -I$(AZTEKAS_PATH)/include -I$(AZTEKAS_PATH)/include/aztekas
else ifeq ($(UNAME_S),Linux)
  # Estamos en Linux
  CXX      = g++
  CXXFLAGS = -Wall -std=c++17 -I$(AZTEKAS_PATH)/include -I$(AZTEKAS_PATH)/include/aztekas
else
  # Otras plataformas, por si acaso quieres un fallback
  CXX      = g++
  CXXFLAGS = -Wall -std=c++17 -I$(AZTEKAS_PATH)/include -I$(AZTEKAS_PATH)/include/aztekas
endif

.PHONY: all clean

all: aztekas-simulation $(EXECUTABLE)

$(EXECUTABLE): $(LOCAL_OBJS) \
               $(COMMON_OBJS) \
               $(GRID_OBJS) \
               $(OBJDIR)/main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(AZTEKAS_PATH)/src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
