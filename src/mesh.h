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

class Vertex;
class Triangle;

class Vertex {
public:
    vector<Triangle*> T;
    double x, y;
    int id;
    int label;
    Vertex(double x, double y, int label, int id);
};

class Triangle {
private:
    double calculate_area();
public:
    vector<Vertex*> V;
    int label;
    int id;
    double area;
    Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id);

};

class BoundEdge {
public:
    vector<Vertex*> V;
    int id;
    int label;
    double length;
    BoundEdge(Vertex* a, Vertex* b, int label, int id);
};

class Mesh {
private:
    friend istream& operator >>(istream &is, Mesh &obj);
    
    
public:
    int Nv;
    int Nt;
    int Ne;
    vector<Vertex> V;
    vector<Triangle> T;
    vector<BoundEdge> E;
};

typedef Mesh maillage;

#endif // __MESH_H__