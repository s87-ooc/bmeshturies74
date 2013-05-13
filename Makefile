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

# -----------------------------------------------------------------------------

# programs solving the assignment

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

# -----------------------------------------------------------------------------
	
# data structures and methods for solving the problem and visualization	

src/algebra.o: src/algebra.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/algebra.o src/algebra.cpp

src/mesh.o: src/mesh.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/mesh.o src/mesh.cpp

src/visualization.o: src/visualization.cpp src/types.h
	g++ $(CXXFLAGS) -c -o src/visualization.o src/visualization.cpp

# -----------------------------------------------------------------------------
	
# unittest

bin/unittest: src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp
	g++ $(CXXFLAGS) -o bin/unittest src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp

test: bin/unittest
	bin/unittest

# -----------------------------------------------------------------------------

# documentation

documentation:
	doxygen docs/_doxygen/projet.doxygen

# -----------------------------------------------------------------------------

# Cleanup targets
	
clean:
	rm -Rf src/*.o bin/*
	
cleandata:
	rm -Rf data/linsys/test_* data/linsys/times_* data/linsys/errors_*
	rm -Rf data/plots/*.p data/plots/*.pdat data/plots/*.mesh data/plots/*.bb data/plots/*.png data/plots/*.avi
	
cleandocs:
	rm -Rf docs/html
	
cleanall: clean cleandata cleandocs