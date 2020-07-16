CC = mpic++
#CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address -std=c++17
CFLAGS = -O3 -std=c++17

OBJECTS = main.o mesh2d.o hydro_flux.o hydro_timestep.o hydro_MH1D.o update.o

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

tonberry: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o tonberry -lm $(LIB)

clean:
	-rm *.o
	-rm tonberry


install: tonberry
	@echo Installing to $(INSTALL_DIR)
	cp tonberry $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/tonberry

main.o: mesh2d.o hydro_flux.o hydro_timestep.o hydro_MH1D.o update.o

update.o: hydro_flux.o hydro_timestep.o hydro_MH1D.o


%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

