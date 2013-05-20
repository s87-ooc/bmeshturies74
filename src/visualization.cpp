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
#include <sstream>
#include <cstdlib>
#include <cstring>

#include <assert.h>
#include <math.h>

#include "algebra.h"
#include "mesh.h"
#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

Plot::Plot() :
mXPtr(0),
mFuncPtr(0)
{
}

Plot::~Plot()
{
}

void Plot::addYVector(Vector& y, const char* args)
{
	mYPtrs.push_back(&y);

	if (args)
	{
		mArgs.push_back(args);
	}
	else
	{
		mArgs.push_back("");
	}
}

void Plot::addScriptLine(const char* line)
{
	if (line)
	{
		mScriptLines.push_back(line);
	}
}

void Plot::setAxisLabel(EPlotAxis axis, const char* label)
{
	if (label)
	{
		if (axis == ePA_X)
		{
			mLabelX = label;
		}
		else
		{
			mLabelY = label;
		}
	}
}

Plot::Plot(const char* name, Vector& x, Vector& y, const char* title, const char* templ, const char* args) :
mName(name),
mXPtr(&x),
mDataType(ePD_VECTOR)
{
	if (title)
	{
		mTitle = title;
	}
	
	if (templ)
	{
		mTemplate = templ;
	}
	
	if (args)
	{
		mArgs.push_back(args);
	}
	else
	{
		mArgs.push_back("");
	}

	mYPtrs.push_back(&y);
}
	
Plot::Plot(const char* name, Vector& x, double (&func)(double), const char* title, const char* templ, const char* args) :
mName(name),
mXPtr(&x),
mFuncPtr(&func),
mDataType(ePD_FUNCTION)
{
	if (title)
	{
		mTitle = title;
	}
	
	if (templ)
	{
		mTemplate = templ;
	}
	
	if (args)
	{
		mArgs.push_back(args);
	}
	else
	{
		mArgs.push_back("");
	}
}

void Plot::generate(EPlotType type, bool run, bool savePNG)
{
	assert(mXPtr);
	
	string fileName = "data/plots/" + mName;
	
	if (type == ePT_GNUPLOT)
	{
		// data file
		{
			ofstream fout((fileName + ".pdat").c_str());
			
			uint dim = mXPtr->size();
			
			if (mDataType == ePD_FUNCTION)
			{
				for (uint i = 0; i < dim; i++)
				{
					fout << (*mXPtr)(i) << " " << mFuncPtr((*mXPtr)(i)) << endl;
				}
			}
			else if (mDataType == ePD_VECTOR)
			{
				assert(dim == mYPtrs[0]->size());
				
				for (uint i = 0; i < dim; i++)
				{
					fout << (*mXPtr)(i) << " " << (*mYPtrs[0])(i);
					for (uint iY = 1; iY < mYPtrs.size(); iY++)
					{
						fout << " " << (*mYPtrs[iY])(i);
					}
					fout << endl;
				}
			}
		}
		
		// script file
		{
			ofstream fout((fileName + ".p").c_str());
			if (mTemplate.empty())
			{
				mTemplate = "data/_gnuplot/default_2d.ptpl";
			}

			ifstream fin(mTemplate.c_str());
			fout << fin.rdbuf() << endl;

			for (uint iLine = 0; iLine < mScriptLines.size(); iLine++)
			{
				fout << mScriptLines[iLine] << endl;
			}

			if (savePNG)
			{
				fout << "set term pngcairo" << endl;
				fout << "set output \"" << fileName << ".png\"" << endl;
			}
			
			if (!mTitle.empty())
			{
				fout << "set title \"" << mTitle << "\"" << endl;
			}
			
			if (!mLabelX.empty())
			{
				fout << "set xlabel \"" << mLabelX << "\"" << endl;
			}
			
			if (!mLabelY.empty())
			{
				fout << "set ylabel \"" << mLabelY << "\"" << endl;
			}
			
			fout << "plot '" << fileName << ".pdat'" << mArgs[0];

			for (uint iY = 1; iY < mYPtrs.size(); iY++)
			{
				fout << ", '" << fileName << ".pdat' using 1:" << iY + 2 ;
				fout << mArgs[iY];
			}

			fout << endl;
			
			if (savePNG)
			{
				fout << "exit gnuplot" << endl;
			}
		}
		
		if (run)
		{
			string cmd = "gnuplot -persist " + fileName + ".p";
			
			if (savePNG)
			{
				cout << " ### rendering " << fileName + ".png" << endl;
			}
			
			system(cmd.c_str());

			if (!savePNG)
			{
				cout << "\n<<< Press any key to continue >>>";
				cin.get();
			}
		}
	}
}

// ----------------------------------------------------------------------------

PlotMesh::PlotMesh() :
mMeshPtr(0),
mFuncPtr(0),
mVectorPtr(0)
{
}

PlotMesh::~PlotMesh()
{
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, const char* title) :
mName(name),
mMeshPtr(&msh),
mDataType(ePD_MESH)
{
	if (title)
	{
		mTitle = title;
	}
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, Vector& vals, const char* title) :
mName(name),
mMeshPtr(&msh),
mVectorPtr(&vals),
mDataType(ePD_VECTOR)
{
	if (title)
	{
		mTitle = title;
	}
}

PlotMesh::PlotMesh(const char* name, Mesh& msh, double (&func)(const Vertex&), const char* title) :
mName(name),
mMeshPtr(&msh),
mFuncPtr(&func),
mDataType(ePD_FUNCTION)
{
	if (title)
	{
		mTitle = title;
	}
}

void PlotMesh::generate(EPlotType type, bool run, bool savePNG, const char* templ, const char* args, uint grid)
{
	assert(mMeshPtr);
	
	string fileName = "data/plots/" + mName;
	
	if (type == ePT_GNUPLOT || type == ePT_GNUPLOT_SURF)
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
				if (type == ePT_GNUPLOT_SURF)
				{
					double d = 1. / (double)grid;
					for (uint i = 0; i <= grid; i++)
					{
						for (uint j = 0; j <= grid; j++)
						{
							fout << i * d << " " << j * d << " " << mMeshPtr->eval(i * d, j * d, *mVectorPtr) << endl;
						}
						if (i < grid)
						{
							fout << endl;
						}
					}
				}
				else
				{
					for (uint t = 0; t < mMeshPtr->countVertices(); t++)
					{
						fout << mMeshPtr->V[t].x << " " << mMeshPtr->V[t].y << " ";
						fout << (*mVectorPtr)(t);
						fout << endl;
					}
				}
			}
		}
		
		// script file
		{
			string sTemplate;
			
			ofstream fout((fileName + ".p").c_str());
			if (!templ || strcmp(templ, "") == 0)
			{
				sTemplate = "data/_gnuplot/default.ptpl";
			}
			else
			{
				sTemplate = templ;
			}

			ifstream fin(sTemplate.c_str());
			fout << fin.rdbuf() << endl;

			if (savePNG)
			{
				fout << "set term pngcairo" << endl;
				fout << "set output \"" << fileName << ".png\"" << endl;
			}
			
			if (!mTitle.empty())
			{
				fout << "set title \"" << mTitle << "\"" << endl;
			}
			
			fout << "splot '" << fileName << ".pdat'";
			
			if (type == ePT_GNUPLOT_SURF)
			{
				fout << " with lines";
			}
			
			fout << endl;
			
			if (savePNG)
			{
				fout << "exit gnuplot" << endl;
			}
		}

		if (run)
		{
			string cmd = "gnuplot -persist " + fileName + ".p";
			
			if (savePNG)
			{
				cout << " ### rendering " << fileName + ".png" << endl;
			}

			system(cmd.c_str());
			
			if (!savePNG)
			{
				cout << "\n<<< Press any key to continue >>>";
				cin.get();
			}
		}
	}
	else if (type == ePT_MEDIT)
	{
		// mesh data file
		{
			ofstream fout((fileName + ".mesh").c_str());
			
			fout << "MeshVersionFormatted 1" << endl;
			fout << endl;
			fout << "Dimension 2" << endl;

			fout << endl;
			fout << "Vertices" << endl;	
			fout << mMeshPtr->countVertices() << endl;
			
			for (uint v = 0; v < mMeshPtr->countVertices(); v++)
			{
				fout << mMeshPtr->V[v].x << " " << mMeshPtr->V[v].y << " 1" << endl;
			}
			
			fout << endl;
			fout << "Triangles" << endl;
			fout << mMeshPtr->countTriangles() << endl;
			
			for (uint t = 0; t < mMeshPtr->countTriangles(); t++)
			{
				fout << mMeshPtr->T[t](0).id + 1 << " " << mMeshPtr->T[t](1).id + 1 << " " << mMeshPtr->T[t](2).id + 1 << " 1" << endl;
			}
			
			fout << endl;
			fout << "Edges" << endl;
			fout << mMeshPtr->countEdges() << endl;
			
			for (uint e = 0; e < mMeshPtr->countEdges(); e++)
			{
				fout << mMeshPtr->E[e](0).id + 1 << " " << mMeshPtr->E[e](1).id + 1 << " 1" << endl;
			}
			
			fout << endl;
			fout << "End" << endl;
			fout << endl;
		}

		// solution file if necessary
		if (mDataType == ePD_FUNCTION || mDataType == ePD_VECTOR)
		{
			ofstream fout((fileName + ".bb").c_str());
			
			fout << "2 1 " << mMeshPtr->countVertices() << " 2" << endl;

			if (mDataType == ePD_FUNCTION)
			{
				for (uint t = 0; t < mMeshPtr->countTriangles(); t++)
				{
					fout << mFuncPtr(mMeshPtr->V[t]) << " ";
				}
			}
			else
			{
				for (uint t = 0; t < mMeshPtr->countTriangles(); t++)
				{
					fout << (*mVectorPtr)(t) << " ";
				}
			}
			
			fout << endl;
		}
		
		if (run)
		{
			string cmd = "medit " + fileName;
			system(cmd.c_str());
		}
	}
}

// ----------------------------------------------------------------------------

void PlotVideo::renderVideo(const char* file, const char* framePrefix, double framerate, const char* title,
		uint width, uint height, EPlotVideo codec)
{
	string cmd = "mencoder mf://";
	cmd += framePrefix;
	cmd += "*.png -mf w=";
	{
		stringstream buf;
		buf << width;
		string sbuf;
		buf >> sbuf;
		cmd += sbuf;
	}
	cmd += ":h=";
	{
		stringstream buf;
		buf << height;
		string sbuf;
		buf >> sbuf;
		cmd += sbuf;
	}
	cmd += ":fps=";
	{
		stringstream buf;
		buf << framerate;
		string sbuf;
		buf >> sbuf;
		cmd += sbuf;
	}
	cmd += ":type=png -ovc copy -oac copy -o ";
	cmd += file;
	
	cout << endl << "### Change into data/plots/ and run" << endl;
	cout << cmd << endl;
	
	//system(cmd.c_str());
}
