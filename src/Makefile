# src/Makefile

# Detectar sistema operativo
OS := $(shell uname -s)

# Opciones por defecto (válidas para ambos sistemas, si no las sobrescribes)
CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Configuraciones específicas para macOS vs. Linux
ifeq ($(OS),Darwin)  # Darwin = macOS
    CXX      = clang++
    CXXFLAGS += -fopenmp -I../include -I../include/aztekas -I.
    # Agrega aquí cualquier otra bandera específica de macOS
    # Por ejemplo, -framework CoreFoundation (si fuera necesario) 
endif

ifeq ($(OS),Linux) 
    CXX      = g++
    CXXFLAGS += -fopenmp -I../include -I../include/aztekas -I.
    # Agrega aquí cualquier otra bandera específica de Linux
    # Por ejemplo, -pthread para threads, -ldl, etc.
endif

TARGET   = main

.PHONY: all clean

all: $(TARGET)

$(TARGET): main.o 
	$(CXX) $(CXXFLAGS) main.o -o $(TARGET)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp 

clean:
	rm -f *.o $(TARGET)