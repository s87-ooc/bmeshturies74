/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: visualization.cpp

 Description: implementation des classes pour la visualization

**************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <assert.h>
#include <math.h>

#include "mesh.h"
#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

Plot::Plot() :
mMesh(0)
{
}

Plot::~Plot()
{
}

Plot::Plot(const char* name, const char* temp, Mesh* msh, EPlotType type, const char* args) :
mName(name),
mTemplate(temp),
mMesh(msh),
mType(type)
{
	if (args)
	{
		mArgs = args;
	}
}

void Plot::generate(bool run)
{
	assert(mMesh);
	
	string fileName = "data/plots/" + mName;
	
	if (mType == ePT_GNUPLOT_SURFACE)
	{
		// data file
		{
			ofstream fout((fileName + ".pdat").c_str());
			
			TVertices::iterator it = mMesh->V.begin();
			
			while (it != mMesh->V.end())
			{
				fout << it->x << " " << it->y << " ";
				
				// TODO: write actual data that's associated to the node here
				fout << pow(it->x, 2) + it->y;
				
				fout << endl;
				
				++it;
			}
		}
		
		// script file
		{
			ofstream fout((fileName + ".p").c_str());
			if (!mTemplate.empty())
			{
				ifstream fin(mTemplate.c_str());
				fout << fin.rdbuf() << endl;
			}
			fout << "splot '" << fileName << ".pdat'";
		}
		
		if (run)
		{
			string cmd = "gnuplot -persist " + fileName + ".p";
			system(cmd.c_str());
		}
	}
	else
	{
		// TODO: define/implement other plot types
	}
}