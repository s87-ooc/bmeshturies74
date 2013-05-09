﻿/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.cpp

 Description: implementation de la classe mesh

**************************************************************/

#include <iostream>
#include <math.h>
#include <assert.h>

#include <map>

#include "mesh.h"


// ----------------------------------------------------------------------------

Vertex::Vertex() :
x(0.),
y(0.),
label(0),
id(0),
T(0)
{
}

Vertex::Vertex(double x, double y, uint label, uint id): 
    x(x), y(y), label(label), id(id), T(0)
{
}

Vertex::Vertex(const Vertex& v)
{
	x = v.x;
	y = v.y;
	label = v.label;
	id = v.id;
    T = v.T;
}

// ----------------------------------------------------------------------------

Triangle::Triangle(): label(-1), id(-1), V(0)
{
    // don't use this constructor
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


BoundEdge::BoundEdge():
V(0), label(-1), id(-1)
{
    // don't use this constructor
}

BoundEdge::BoundEdge(Vertex* a, Vertex* b, int label, int id):
label(label), id(id), V(2)
{
    V[0]=a;
    V[1]=b;
    length = calculate_length();
    edgeOf = findTriangle();
}

double BoundEdge::calculate_length()
{
    return sqrt( pow(V[0]->x - V[1]->x, 2) + pow(V[0]->y - V[1]->y, 2));
}

Triangle* BoundEdge::findTriangle() const
{

    for( uint k = 0; k < V[0]->T.size(); k++)
    {
        for( uint l = 0; l < V[1]->T.size(); l++)
        {
            if( V[0]->T[k]->id == V[1]->T[l]->id)
            {
                cout << V[0]->T[k]->id << endl;
                return V[0]->T[k];
            }
        }
    }
    assert(false);
}

bool BoundEdge::inEdge(const Vertex* v) const
{
   return ( v->id == V[0]->id || v->id == V[1]->id );
}

Vertex& BoundEdge::findOppositeVertex() const
{
    //Triangle* edgeOf = findTriangle();
    cout << "looking for opposite vertex in triangle " << edgeOf->id << endl;
    for(int i=0; i<3; i++)
    {
        if ( !inEdge(edgeOf->V[i]) ) 
        {
            return *(edgeOf->V[i]);
        }
    }
    assert(false);
}

Vertex& BoundEdge::operator() (uint vertex) const
{
    assert( vertex >=0 && vertex <=2);
    return *V[vertex];
}

istream& operator>>(istream &is, Mesh &M)
{
	// // we're resetting the mesh, clear stuff if needed
	if (M.V)
	{
		delete[] M.V;
	}

    // lecture du nombre de sommets, aretes et triangles
    is >> M.Nv >> M.Nt >> M.Ne;
    
    M.V = new Vertex[M.Nv];
    M.T = new Triangle[M.Nt];
    M.E = new BoundEdge[M.Ne];
	
    // initialiser tableau pour stocker les sommets, aretes et triangles
    
    for(int k = 0; k < M.Nv; k++)
    {
        double x, y; int label;
        is >> x >> y >> label;

        // the ids are 0 indexed
        
        Vertex v(x, y, label, k);
		//Vertex v(x + 0.0001, y + 0.0001, label, k);
		M.V[k] = v;
    }

    for(int k = 0; k < M.Nt; k++)
    {
        int label;
        int a, b, c;
        is >> a >> b >> c >> label;

        Triangle t(&M.V[a-1], &M.V[b-1], &M.V[c-1], label, k);

        M.T[k] = t;

        M.V[a-1].T.push_back(&M.T[k]);
        M.V[b-1].T.push_back(&M.T[k]);
        M.V[c-1].T.push_back(&M.T[k]);

    }

    for(int k = 0; k < M.Ne; k++)
    {
        int label;
        int a, b;
        is >> a >> b >> label;

        BoundEdge e(&M.V[a-1], &M.V[b-1], label, k);
        M.E[k] = e;
    }


    // for(int k = 0; k < M.Nv; k++)
    // {
    //     cout << "Vertex " << M.V[k].id << " is in triangles";
    //     for(int l = 0; l < M.V[k].T.size(); l++)
    //     {
    //         cout << " " << M.V[k].T[l]->id;
    //     }
    //     cout << endl;
    // }
    
    cout << "Finished reading Mesh from file" << endl;
    
    return is;
}

Mesh::Mesh(const char* filename)
{
    ifstream meshfile;
    meshfile.open(filename, ifstream::in);
    
    meshfile >> (*this);
}

Mesh::~Mesh()
{
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