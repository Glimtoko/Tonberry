CC = mpic++
#CFLAGS = -g -Wall -Wextra -ansi -pedantic -fsanitize=address -std=c++17
CFLAGS = -O3 -std=c++17

OBJECTS = main.o mesh2d.o flux.o update.o

BOOST_DIR=/prod/boost
BOOST_INC=$(BOOST_DIR)/include
BOOST_LIB=$(BOOST_DIR)/lib

INC = -I$(BOOST_INC)
LIB = -L$(BOOST_LIB) -lboost_program_options -lboost_mpi -lnetcdf_c++4 -lnetcdf

.PHONY: clean

tonberry: $(OBJECTS)
	$(CC) $(CFLAGS) *.o -o tonberry -lm $(LIB)

clean:
	-rm *.o
	-rm tonberry


main.o: mesh2d.o flux.o update.o

update.o: flux.o



%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

