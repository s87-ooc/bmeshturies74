/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.cpp

 Description: implementation de la classe mesh

**************************************************************/

#include <iostream>
#include <math.h>

#include <map>

#include "mesh.h"

using namespace std;

// ----------------------------------------------------------------------------

Vertex::Vertex(double x, double y, int label, int id): 
    x(x), y(y), label(label), id(id)
{
}

Triangle::Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id):
    label(label), id(id)
{
    V.push_back(a);
    V.push_back(b);
    V.push_back(c);
    area = calculate_area();
}

double Triangle::calculate_area()
{  
    double a,b;
    a = (V[0]->x - V[2]->x)*(V[1]->y - V[0]->y);
    b = (V[0]->x - V[1]->x)*(V[2]->y - V[0]->y);
    return 0.5*fabs(a-b);
}

Vertex& Triangle::operator() (uint vertex) const
{
    assert( vertex >=0 && vertex <=2);
    return *V[vertex];
}


BoundEdge::BoundEdge(Vertex* a, Vertex* b, int label, int id):
label(label), id(id)
{
    V.push_back(a);
    V.push_back(b);
    length = calculate_length();
}

double BoundEdge::calculate_length()
{
    return sqrt( pow(V[0]->x - V[1]->x, 2) + pow(V[0]->y - V[1]->y, 2));
}

Vertex& BoundEdge::operator() (uint vertex) const
{
    assert( vertex >=0 && vertex <=2);
    return *V[vertex];
}

istream& operator>>(istream &is, Mesh &M)
{
    // lecture du nombre de sommets, aretes et triangles
    is >> M.Nv >> M.Nt >> M.Ne;
    
    
    // initialiser tableau pour stocker les sommets, aretes et triangles
    
    for(int k = 0; k < M.Nv; k++)
    {
        double x, y; int label;
        is >> x >> y >> label;

        // the ids are 0 indexed
        Vertex v(x, y, label, k);
        M.V.push_back(v);
    }

    for(int k = 0; k < M.Nt; k++)
    {
        int label;
        int a, b, c;
        is >> a >> b >> c >> label;

        Triangle t(&M.V[a-1], &M.V[b-1], &M.V[c-1], label, k);

        M.T.push_back(t);

        for(int i=0; i<3; i++)
        {
            M.T[k].V[i]->T.push_back(&(M.T[k]));
        }

    }
    for(int k = 0; k < M.Ne; k++)
    {
        int label;
        int a, b;
        is >> a >> b >> label;

        BoundEdge e(&(M.V[a-1]), &(M.V[b-1]), label, k);
        M.E.push_back(e);
    }
    
    return is;
}

Mesh::Mesh(const char* filename)
{
    ifstream meshfile;
    meshfile.open(filename, ifstream::in);
    
    meshfile >> (*this);
}

uint Mesh::countVertices() const
{
    return Nv;
}
uint Mesh::countTriangles() const
{
    return Nt;
}
uint Mesh::countEdges() const
{
    return Ne;
}