###############################################################################
#                                                                             #
# Projet                                                                      #
#                                                                             #
# (C) 2013 Charles Podkanski (charles@gmail.com),                             #
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

all: bin/helmholtz bin/wave

# solutions of the two parts
bin/helmholtz: src/algebra.o src/mesh.o src/visualization.o src/helmholtz.cpp
	g++ $(CXXFLAGS) -o bin/helmholtz src/algebra.o src/mesh.o src/visualization.o src/helmholtz.cpp

bin/wave: src/algebra.o src/mesh.o src/visualization.o src/wave.cpp
	g++ $(CXXFLAGS) -o bin/wave src/algebra.o src/mesh.o src/visualization.o src/wave.cpp

# data structures and methods for solving the problem and visualization	
src/algebra.o: src/algebra.cpp
	g++ $(CXXFLAGS) -c -o src/algebra.o src/algebra.cpp

src/mesh.o: src/mesh.cpp
	g++ $(CXXFLAGS) -c -o src/mesh.o src/mesh.cpp

src/visualization.o: src/visualization.cpp
	g++ $(CXXFLAGS) -c -o src/visualization.o src/visualization.cpp

# tests the implementation of the classes
bin/unittest: src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp
	g++ $(CXXFLAGS) -o bin/unittest src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp

test: bin/unittest
	bin/unittest
	
clean:
	rm -Rf src/*.o bin/*