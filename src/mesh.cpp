/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.cpp

 Description: implementation de la classe mesh

**************************************************************/

#include <iostream>
#include <math.h>

#include <map>

#include "mesh.h"


// ----------------------------------------------------------------------------

Vertex::Vertex() :
x(0.),
y(0.),
label(0),
id(0)
{
}

Vertex::Vertex(double x, double y, uint label, uint id): 
    x(x), y(y), label(label), id(id)
{
}

Vertex::Vertex(const Vertex& v)
{
	x = v.x;
	y = v.y;
	label = v.label;
	id = v.id;
}

// ----------------------------------------------------------------------------

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
    edgeOf = findTriangle();
}

double BoundEdge::calculate_length()
{
    return sqrt( pow(V[0]->x - V[1]->x, 2) + pow(V[0]->y - V[1]->y, 2));
}

Triangle* BoundEdge::findTriangle()
{
    TTrianglesP& T0 = V[0]->T;
    TTrianglesP& T1 = V[1]->T;

    for( uint k = 0; k < T0.size(); k++)
    {
        for( uint l = 0; l < T1.size(); l++)
        {
            if( T0[k]->id == T1[l]->id)
            {
                return T0[k];
            }
        }
    }
    assert(false);
}

bool BoundEdge::inEdge(const Vertex* v) const
{
    return ( v->id == V[0]->id || v->id == V[1]->id );
}

Vertex& BoundEdge::findOppositeVertex()
{
    for(int i=0; i<3; i++)
    {
        if ( !inEdge(edgeOf->V[i]) ) 
        {
            return *V[i];
        }
    }
    assert(false);
}

// Vector& BoundEdge::normal()
// {
//     Vector n(2);
//     Vertex& opp = findOppositeVertex();

//     // calculate normal to edge
//     n(0) = V[1]->x - V[0]->x;
//     n(1) = V[1]->y - V[0]->y;

//     // check orientation
//     double D = n(0)*( opp.x - V[0]->x ) + n(1)*( opp.y - V[0]->y);

//     if (D > 0)
//     {
//         // UGLY, could use - operator ?
//         n(0) = -n(0);
//         n(1) = -n(1);
//     }

//     // normalize
//     n *= 1./n.norm2();

//     return n;
// }


Vertex& BoundEdge::operator() (uint vertex) const
{
    assert( vertex >=0 && vertex <=2);
    return *V[vertex];
}

istream& operator>>(istream &is, Mesh &M)
{
	// we're resetting the mesh, clear stuff if needed
	if (M.V)
	{
		delete[] M.V;
	}

    // lecture du nombre de sommets, aretes et triangles
    is >> M.Nv >> M.Nt >> M.Ne;
    
	// allocate memory
	M.V = new Vertex[M.Nv];
	
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

Mesh::Mesh(const char* filename) :
V(0)
{
    ifstream meshfile;
    meshfile.open(filename, ifstream::in);
    
    meshfile >> (*this);
}

Mesh::~Mesh()
{
	if (V)
	{
		delete[] V;
	}
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