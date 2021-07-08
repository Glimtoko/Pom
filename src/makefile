CC = g++
#CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address -std=c++17
#CFLAGS = -g
CFLAGS = -O3 -std=c++17 -fopenmp

OBJECTS = main.o mesh2d.o hydro_flux.o hydro_timestep.o  hydro_MH2D.o

# Installation directory
INSTALL_DIR := /prod/Tonberry/bin/

# Expected library locations
PROD := /prod
SILO_DIR := $(PROD)/silo
TIO_DIR := $(PROD)/typhonio
HDF5_DIR := $(PROD)/HDF5_1.12.0_ser_gnu

INC = -I$(SILO_DIR)/include -I$(TIO_DIR)/include -I$(HDF5_DIR)/include
LIB = -L$(SILO_DIR)/lib -lsilo -L$(TIO_DIR)/lib -ltyphonio -L$(HDF5_DIR)/lib -lhdf5 -Wl,-rpath $(HDF5_DIR)/lib

VPATH=.:./hydro:./mesh

.PHONY: clean

pom: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o pom -lm $(LIB)

clean:
	-rm *.o
	-rm pom


install: pom
	@echo Installing to $(INSTALL_DIR)
	cp pom $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/pom

main.o: mesh2d.o hydro_flux.o hydro_timestep.o hydro_MH2D.o


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

