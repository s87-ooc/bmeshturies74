###############################################################################
#                                                                             #
# Projet                                                                      #
#                                                                             #
# (C) 2013 Charles Podkanski (charles@podkanski.com),                         #
#          Stjepan Stamenkovic (stjepan@stjepan.net)                          #
#                                                                             #
# ---                                                                         #
#                                                                             #
#     Fichier: Makefile                                                       #
#                                                                             #
# Description: Contrôle la compilation                                        #
#                                                                             #
###############################################################################

# TODO: Add CXX include paths, etc
# TODO: Debug / Optimize builds

# DEFAULT
CXXFLAGS=

# DEBUG
#CXXFLAGS=-g

# OPTIMIZE
#CXXFLAGS=-O3

all: bin/helmholtz bin/wave bin/solvertest bin/unittest

# solutions of the two parts
helmholtz: bin/helmholtz
	bin/helmholtz
	
wave: bin/wave
	bin/wave

solvertest: bin/solvertest
	bin/solvertest
	
bin/helmholtz: src/algebra.o src/mesh.o src/visualization.o src/helmholtz.cpp
	g++ $(CXXFLAGS) -o bin/helmholtz src/algebra.o src/mesh.o src/visualization.o src/helmholtz.cpp

bin/wave: src/algebra.o src/mesh.o src/visualization.o src/wave.cpp
	g++ $(CXXFLAGS) -o bin/wave src/algebra.o src/mesh.o src/visualization.o src/wave.cpp
	
bin/solvertest: src/algebra.o src/mesh.o src/visualization.o src/solvertest.cpp
	g++ $(CXXFLAGS) -o bin/solvertest src/algebra.o src/mesh.o src/visualization.o src/solvertest.cpp

# data structures and methods for solving the problem and visualization	
src/algebra.o: src/algebra.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/algebra.o src/algebra.cpp

src/mesh.o: src/mesh.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/mesh.o src/mesh.cpp

src/visualization.o: src/visualization.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/visualization.o src/visualization.cpp

# tests the implementation of the classes
bin/unittest: src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp
	g++ $(CXXFLAGS) -o bin/unittest src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp

test: bin/unittest
	bin/unittest
	
clean:
	rm -Rf src/*.o bin/*