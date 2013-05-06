/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
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

#include "algebra.h"
#include "mesh.h"
#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

PlotMesh::PlotMesh() :
mMeshPtr(0)
{
}

PlotMesh::~PlotMesh()
{
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, EPlotType type, const char* templ, const char* args) :
mName(name),
mMeshPtr(&msh),
mDataType(ePD_MESH),
mType(type)
{
	if (templ)
	{
		mTemplate = templ;
	}
	
	if (args)
	{
		mArgs = args;
	}
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, Vector& vals, EPlotType type, const char* templ, const char* args) :
mName(name),
mMeshPtr(&msh),
mVectorPtr(&vals),
mDataType(ePD_VECTOR),
mType(type)
{
	if (templ)
	{
		mTemplate = templ;
	}
	
	if (args)
	{
		mArgs = args;
	}
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, double (&func)(const Vertex&), EPlotType type, const char* templ, const char* args) :
mName(name),
mMeshPtr(&msh),
mFuncPtr(&func),
mDataType(ePD_FUNCTION),
mType(type)
{
	if (templ)
	{
		mTemplate = templ;
	}
	
	if (args)
	{
		mArgs = args;
	}
}

void PlotMesh::generate(bool run)
{
	assert(mMeshPtr);
	
	string fileName = "data/plots/" + mName;
	
	if (mType == ePT_GNUPLOT)
	{
		// data file
		{
			ofstream fout((fileName + ".pdat").c_str());
			
			//TVertices::iterator it = mMeshPtr->V.begin();
			
			if (mDataType == ePD_MESH)
			{
				for (uint t = 0; t < mMeshPtr->countVertices(); t++)
				{
					fout << mMeshPtr->V[t].x << " " << mMeshPtr->V[t].y << " ";
					fout << endl;
				}
			}
			else if (mDataType == ePD_FUNCTION)
			{
				for (uint t = 0; t < mMeshPtr->countVertices(); t++)
				{
					fout << mMeshPtr->V[t].x << " " << mMeshPtr->V[t].y << " ";
					fout << mFuncPtr(mMeshPtr->V[t]);
					fout << endl;
				}
			}
			else if (mDataType == ePD_VECTOR)
			{
				for (uint t = 0; t < mMeshPtr->countVertices(); t++)
				{
					fout << mMeshPtr->V[t].x << " " << mMeshPtr->V[t].y << " ";
					fout << (*mVectorPtr)(t);
					fout << endl;
				}
			}
		}
		
		// script file
		{
			ofstream fout((fileName + ".p").c_str());
			if (mTemplate.empty())
			{
				mTemplate = "data/_gnuplot/default.ptpl";
			}

			ifstream fin(mTemplate.c_str());
			fout << fin.rdbuf() << endl;

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