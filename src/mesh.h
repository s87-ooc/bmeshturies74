#ifndef __MESH_H__
#define __MESH_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.h

 Description: déclaration de la classe mesh

**************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "types.h"

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
private:
	/*
	vector<Triangle*> mT;
	double mX, mY;
    unsigned int mId;
    unsigned int mLabel;
	*/
public:
    TTrianglesP T;
    double x, y;
    int id; // redundant
    int label;
    Vertex(double x, double y, int label, int id);
};

// ----------------------------------------------------------------------------

class Triangle {
private:
    double calculate_area();
	/*
	vector<Vertex*> mV;
    int mLabel;
    int mId;
    double mArea;
	*/
public:
	// Vertex& V[3];
    TVerticesP V;
    int label;
    int id;
    double area;
    Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id);

};

// ----------------------------------------------------------------------------

class BoundEdge {
private:
	/*
	int mId;
    int mLabel;
    double mLength;
	*/
public:
	// Vertex& V[2];
    TVerticesP V;
    int id;
    int label;
    double length;
    BoundEdge(Vertex* a, Vertex* b, int label, int id);
};

// ----------------------------------------------------------------------------

class Plot;

class Mesh {
private:
    friend std::istream& operator >>(std::istream &is, Mesh &obj);
    friend class Plot;
	
public:
	// Stjepan: the following ints should be only local in >> and UNSIGNED!!! :-P
    // unsigned int countVertices() const;
	// unsigned int countTriangles() const;
	// unsigned int countEdges() const;
	int Nv;
    int Nt;
    int Ne;
	
    TVertices V;
    TTriangles T;
    TEdges E;
};

#endif // __MESH_H__