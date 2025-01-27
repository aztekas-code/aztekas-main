# src/Makefile
include $(AZTEKAS_PATH)/src/common/Makefile

# Si no se define desde afuera, será "."
OBJDIR ?= .
CXX      = g++
CXXFLAGS = -Wall -std=c++17 -I$(AZTEKAS_PATH)/include/aztekas

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(LOCAL_OBJS) \
               $(COMMON_OBJS) \
               $(OBJDIR)/main.o 
	$(CXX) $(CXXFLAGS) -o $@ $^ && mv *.o $(OBJDIR)/obj

$(OBJDIR)/%.o: $(AZTEKAS_PATH)/src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
