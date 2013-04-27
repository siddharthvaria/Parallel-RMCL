CC=mpicxx
CFLAGS=-c -Wall -g -O0

all: rmcl

rmcl: spMat_utilities.o spMat_multiply.o rmcl_utilities.o rmcl.o
	$(CC) spMat_utilities.o spMat_multiply.o rmcl_utilities.o rmcl.o -g -O0 -o rmcl_mpi

spMat_utilities.o: spMat_utilities.cpp spMat.hpp rmcl.hpp
	$(CC) $(CFLAGS) spMat_utilities.cpp

spMat_multiply.o: spMat_multiply.cpp spMat_utilities.o
	$(CC) $(CFLAGS) spMat_multiply.cpp
	
rmcl_utilities.o: rmcl_utilities.cpp spMat_multiply.o
	$(CC) $(CFLAGS) rmcl_utilities.cpp

rmcl.o: rmcl.cpp rmcl_utilities.o
	$(CC) $(CFLAGS) rmcl.cpp

clean:
	rm -rf *o rmcl_mpi
