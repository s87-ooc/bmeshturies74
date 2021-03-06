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
#include "algebra.h"


// ----------------------------------------------------------------------------

Vertex::Vertex() :
x(0.), y(0.), label(0), id(0), T(0)
{
}

Vertex::Vertex(double x, double y, uint label, uint id): 
x(x), y(y), label(label), id(id), T(0)
{
}

Vertex::Vertex(const Vertex& v) :
x(v.x), y(v.y), label(v.label), id(v.id), T(v.T)
{
}

ostream& operator<<(ostream& os, const Vertex& v)
{
    os << "Vertex: ID:" << v.id << " x:" << v.x << " y:" << v.y << endl;
    return os;
}

double distance(Vertex* v, Vertex* w)
{
    return sqrt( pow(v->x - w->x,2) + pow(v->y - w->y, 2));
}

// ----------------------------------------------------------------------------

Triangle::Triangle():
label(-1),
id(-1),
V(0)
{
    // don't use this constructor
}

Triangle::Triangle(Vertex* a, Vertex* b, Vertex* c, int label, int id):
label(label),
id(id)
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

double Triangle::circumcircleDiameter() const
{
    double A = distance(V[0], V[1]);
    double B = distance(V[1], V[2]);
    double C = distance(V[2], V[0]);

    return A*B*C / (2.0 * area);
}

double Triangle::incircleDiameter() const
{
    double A = distance(V[0], V[1]);
    double B = distance(V[1], V[2]);
    double C = distance(V[2], V[0]);

    return 4.0 * area / (A+B+C);
}

Vertex& Triangle::operator() (uint vertex) const
{
    assert(vertex >=0 && vertex <=2);
    return *V[vertex];
}

ostream& operator<<(ostream& os, const Triangle& t)
{
    os  << "Triangle: ID:" << t.id << endl
        << "|-" << *t.V[0]
        << "|-" << *t.V[1]
        << "|-" << *t.V[2];
    return os;
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

ostream& operator<<(ostream& os, const BoundEdge& e)
{
    os  << "BoundEdge: ID:" << e.id << endl
        << "|-" << *e.V[0]
        << "|-" << *e.V[1];
    return os;
}

istream& operator>>(istream &is, Mesh &M)
{
	// we're resetting the mesh, clear stuff if needed
	SAFE_ARRDELETE(M.V);
	SAFE_ARRDELETE(M.T);
	SAFE_ARRDELETE(M.E);

    // lecture du nombre de sommets, aretes et triangles
    is >> M.Nv >> M.Nt >> M.Ne;

	// initialiser tableau pour stocker les sommets, aretes et triangles
    M.V = new Vertex[M.Nv];
    M.T = new Triangle[M.Nt];
    M.E = new BoundEdge[M.Ne];	
	
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
    
//    cout << "Finished reading Mesh from file" << endl;
    
    return is;
}

Mesh::Mesh() :
V(0),
T(0),
E(0),
Nv(0),
Nt(0),
Ne(0)
{
}
/*
Mesh::Mesh(const Mesh& mesh)
{
	SAFE_ARRDELETE(V);
	SAFE_ARRDELETE(T);
	SAFE_ARRDELETE(E);

	Nv = mesh.countVertices();
	Nt = mesh.countTriangles();
	Ne = mesh.countEdges();

	V = new Vertex[mesh.countVertices()];
	T = new Triangle[mesh.countTriangles()];
	E = new BoundEdge[mesh.countEdges()];

	for (uint i = 0; i < mesh.countVertices(); i++)
	{
		V[i] = mesh.V[i];
	}

	for (uint i = 0; i < mesh.countTriangles(); i++)
	{
		T[i] = mesh.T[i];
	}

	for (uint i = 0; i < mesh.countEdges(); i++)
	{
		E[i] = mesh.E[i];
	}
}

Mesh& Mesh::operator=(const Mesh& mesh)
{
	SAFE_ARRDELETE(V);
	SAFE_ARRDELETE(T);
	SAFE_ARRDELETE(E);

	Nv = mesh.countVertices();
	Nt = mesh.countTriangles();
	Ne = mesh.countEdges();

	V = new Vertex[mesh.countVertices()];
	T = new Triangle[mesh.countTriangles()];
	E = new BoundEdge[mesh.countEdges()];

	for (uint i = 0; i < mesh.countVertices(); i++)
	{
		V[i] = mesh.V[i];
	}

	for (uint i = 0; i < mesh.countTriangles(); i++)
	{
		T[i] = mesh.T[i];
	}

	for (uint i = 0; i < mesh.countEdges(); i++)
	{
		E[i] = mesh.E[i];
	}
}*/

Mesh::~Mesh()
{
	SAFE_ARRDELETE(V);
	SAFE_ARRDELETE(T);
	SAFE_ARRDELETE(E);
}

Mesh::Mesh(const char* filename) :
V(0),
T(0),
E(0)
{
    ifstream meshfile;
    meshfile.open(filename, ifstream::in);
    
    meshfile >> (*this);
}

Mesh::Mesh(uint Nv, uint Nt, uint Ne):
Nv(Nv), Nt(Nt), Ne(Ne)
{
    V = new Vertex[Nv];
    T = new Triangle[Nt];
    E = new BoundEdge[Ne];
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

// TODO: implement this more elegant

#define DIST(u,v) pow(u.x - v.x, 2) + pow(u.y - v.y, 2)

double Mesh::maxCircumcircleDiameter() const
{
	double maxDiam = 0.0;
	
	for (uint i = 0; i < Nt; i++)
	{
		maxDiam = max(maxDiam, T[i].circumcircleDiameter() );
	}
	
	return maxDiam;
}

double Mesh::maxIncircleDiameter() const
{
    double maxDiam = 0.0;
    
    for (uint i = 0; i < Nt; i++)
    {
        maxDiam = max(maxDiam, T[i].incircleDiameter() );
    }
    
    return maxDiam;
}

double Mesh::eval(double x, double y, const Vector& uh)
{
	assert(uh.size() == countVertices());
	
	Triangle* TPtr = 0;
	
	double c0, c1, c2;
		
	// find the triangle the point is in, then coordinates are in [0, 1]
	// we're tolerating deviations of EQ_TOL to evaluate especially corners of triangles correctly
	for (uint t = 0; t < countTriangles(); t++)
	{
		c0 = ((T[t](1).y - T[t](2).y) * (x - T[t](2).x) + (T[t](2).x - T[t](1).x) * (y - T[t](2).y)) / (2 * T[t].area);
		c1 = ((T[t](2).y - T[t](0).y) * (x - T[t](2).x) + (T[t](0).x - T[t](2).x) * (y - T[t](2).y)) / (2 * T[t].area);
		if (c0 <= 1 + EQ_TOL && c0 >= 0 - EQ_TOL && c1 <= 1 + EQ_TOL && c1 >= 0 - EQ_TOL)
		{
			c2 = 1 - c0 - c1;
			if (c2 <= 1 + EQ_TOL && c2 >= 0 - EQ_TOL)
			{
				TPtr = &T[t];
				break;
			}
		}
	}
	
	if (!TPtr)
	{
		cout << "Point (" << x << ", " << y << ") is not in the mesh!" << endl;
		return 0.;
	}
	
	return c0 * uh((*TPtr)(0).id) + c1 * uh((*TPtr)(1).id) + c2 * uh((*TPtr)(2).id);
}
