#ifndef __MESH_H__
#define __MESH_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.h

 Description: déclaration de la classe mesh

**************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>

#include "types.h"

using namespace std;

// ----------------------------------------------------------------------------

class Vertex;
class Triangle;
class BoundEdge;

typedef std::vector<Vertex> TVertices;
typedef std::vector<Vertex*> TVerticesP;

typedef std::vector<Triangle> TTriangles;
typedef std::vector<Triangle*> TTrianglesP;

typedef std::vector<BoundEdge> TEdges;

// ----------------------------------------------------------------------------

class Vertex;
class Triangle;

class Vertex {
public:
    TTrianglesP T;
	
    double x, y;
    
	uint id; // redundant
    uint label;
	
	Vertex();
	Vertex(const Vertex& v);    
	Vertex(double x, double y, uint label, uint id);
};

// ----------------------------------------------------------------------------

class Triangle {
private:
    double calculate_area();
public:
	// Vertex& V[3];
    TVerticesP V;
    int label;
    int id;
    double area;
    Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id);
    
    Vertex& operator() (uint vertex) const;
};

// ----------------------------------------------------------------------------

class BoundEdge {
private:
    Triangle* findTriangle();
    double calculate_length();
public:
	// Vertex& V[2];
    TVerticesP V;
    int id;
    int label;
    double length;
    Triangle* edgeOf;
    BoundEdge(Vertex* a, Vertex* b, int label, int id);

    Vertex& operator() (uint vertex) const;

};

// ----------------------------------------------------------------------------

class Plot;

class Mesh {
private:
    uint Nv;
    uint Nt;
    uint Ne;

    friend std::istream& operator >>(std::istream &is, Mesh &obj);
    friend class Plot;
	
public:
    Mesh(const char* filename);
	~Mesh();
	
    uint countVertices() const;
	uint countTriangles() const;
	uint countEdges() const;
	
    //TVertices V;
	Vertex* V;
    TTriangles T;
    TEdges E;
};

#endif // __MESH_H__