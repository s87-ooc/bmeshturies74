/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: mesh.cpp

 Description: implementation de la classe mesh

**************************************************************/

#include <iostream>

#include "mesh.h"

using namespace std;

// ----------------------------------------------------------------------------

CMesh::CMesh()
{
	cout << "Mesh created" << endl;
}

CMesh::~CMesh()
{
	cout << "Mesh destroyed" << endl;
}