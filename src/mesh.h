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

using namespace std;

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
    vector<Triangle*> T;
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
    vector<Vertex*> V;
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
    vector<Vertex*> V;
    int id;
    int label;
    double length;
    BoundEdge(Vertex* a, Vertex* b, int label, int id);
};

// ----------------------------------------------------------------------------

class Plot;

typedef vector<Vertex> TVertices;
typedef vector<Triangle> TTriangles;
typedef vector<BoundEdge> TEdges;

class Mesh {
private:
    friend istream& operator >>(istream &is, Mesh &obj);
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