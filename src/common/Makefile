# src/common/Makefile

OBJDIR ?= .

COMMON_SRCS = messages.cpp \
		boundaries.cpp
COMMON_OBJS = $(COMMON_SRCS:%.cpp=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: $(AZTEKAS_PATH)/src/common/%.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@
