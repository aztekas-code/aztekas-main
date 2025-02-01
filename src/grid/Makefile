# src/common/Makefile

OBJDIR ?= .

GRID_SRCS = allocate.cpp
GRID_OBJS = $(GRID_SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: $(AZTEKAS_PATH)/src/grid/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
