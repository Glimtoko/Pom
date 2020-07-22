CC = g++
CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address -std=c++17
#CFLAGS = -g
#CFLAGS = -O3 -std=c++17 -fopenmp

OBJECTS = main.o mesh2d.o hydro_flux.o hydro_timestep.o hydro_MH1D.o hydro_MH2D.o

# Installation directory
INSTALL_DIR := /prod/Tonberry/bin/

# Expected library locations
PROD := /usr/local/prod
SILO_DIR := $(PROD)/silo-4.10.2-bsd
HDF5_DIR := $(PROD)/hdf5_oldapi

INC = -I$(SILO_DIR)/include -I$(HDF5_DIR)/include
LIB = -L$(SILO_DIR)/lib -lsiloh5 -lnetcdf_c++4 -L$(HDF5_DIR)/lib -lhdf5 -Wl,-rpath $(PROD)/hdf5/lib

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

main.o: mesh2d.o hydro_flux.o hydro_timestep.o hydro_MH1D.o hydro_MH2D.o


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

