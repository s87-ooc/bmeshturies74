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
typedef std::vector<BoundEdge*> TEdgesP;


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

    friend ostream& operator<<(ostream& os, const Vertex& v);

};

double distance(const Vertex* v, const Vertex* w);

// ----------------------------------------------------------------------------

class Triangle {
private:
    double calculate_area();
public:
    double circumcircleDiameter() const;
    double incircleDiameter() const;
	// Vertex& V[3];
    TVerticesP V;
    int label;
    int id;
    double area;
    Triangle();
    Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id);
    
    Vertex& operator() (uint vertex) const;

    friend ostream& operator<<(ostream& os, const Triangle& t);

};

// ----------------------------------------------------------------------------

class BoundEdge {
private:
    Triangle* findTriangle() const;
    double calculate_length();
public:
	// Vertex& V[2];
    TVerticesP V;
    int id;
    int label;
    double length;
    Triangle* edgeOf;
    BoundEdge();
    BoundEdge(Vertex* a, Vertex* b, int label, int id);

    Vertex& findOppositeVertex() const;
    Vertex& operator() (uint vertex) const;
    bool inEdge(const Vertex* v) const;

    friend ostream& operator<<(ostream& os, const BoundEdge& e);


};

// ----------------------------------------------------------------------------

class Plot;
class Vector;

class Mesh {
private:
    uint Nv;
    uint Nt;
    uint Ne;

    friend std::istream& operator >>(std::istream &is, Mesh &obj);
    friend class Plot;
	
public:
    Mesh(const char* filename);
    Mesh(uint Nv, uint Nt, uint Ne);
	~Mesh();
	
    uint countVertices() const;
	uint countTriangles() const;
	uint countEdges() const;
	
	/** return the maximal diameter of a triangle in this mesh */
	double maxCircumcircleDiameter() const;
    double maxIncircleDiameter() const;

	
	/** evaluate a function associated to the discrete space of the mesh with coefficients from uh */
	double eval(double x, double y, const Vector& uh); 
	
	Vertex* V;
    Triangle* T;
    BoundEdge* E;
};

#endif // __MESH_H__