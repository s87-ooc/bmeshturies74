# TODO: Add CXX include paths, etc
# TODO: Debug / Optimize builds

# DEFAULT
CXXFLAGS=

# DEBUG
#CXXFLAGS=-g

# OPTIMIZE
#CXXFLAGS=-O3

all: helmholtz wave

# solutions of the two parts
helmholtz: algebra.o mesh.o visualization.o src/helmholtz.cpp
	g++ $(CXXFLAGS) -o bin/helmholtz src/algebra.o src/mesh.o src/visualization.o src/helmholtz.cpp

wave: algebra.o mesh.o visualization.o src/wave.cpp
	g++ $(CXXFLAGS) -o bin/wave src/algebra.o src/mesh.o src/visualization.o src/wave.cpp

# data structures and methods for solving the problem and visualization	
algebra.o: src/algebra.cpp
	g++ $(CXXFLAGS) -c -o src/algebra.o src/algebra.cpp

mesh.o: src/mesh.cpp
	g++ $(CXXFLAGS) -c -o src/mesh.o src/mesh.cpp

visualization.o: src/visualization.cpp
	g++ $(CXXFLAGS) -c -o src/visualization.o src/visualization.cpp

# tests the implementation of the classes
unittest: algebra.o mesh.o visualization.o src/unittest.cpp
	g++ $(CXXFLAGS) -o bin/unittest src/algebra.o src/mesh.o src/visualization.o src/unittest.cpp

test: unittest
	bin/unittest
	
clean:
	rm -Rf src/*.o bin/*