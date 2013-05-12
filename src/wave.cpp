/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenković (stjepan@stjepan.net)

 ---

     Fichier: wave.cpp

 Description: solution de la deuxième partie

**************************************************************/

#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iomanip>

// ---

#include "algebra.h"
#include "mesh.h"
#include "visualization.h"
#include "clock.h"

using namespace std;

// ----------------------------------------------------------------------------

// @@@
template<class T>
void _writeMatrix(const char* fn, const T& m)
{
	ofstream os(fn);
	for (uint i = 0; i < m.sizeRows(); i++)
	{
		for (uint j = 0; j < m.sizeColumns(); j++)
		{
			os << m(i, j) << endl;
		}
	}
}

template<class T>
void _writeVector(const char* fn, const T& v)
{
	ofstream os(fn);
	for (uint i = 0; i < v.size(); i++)
	{
		os << v(i) << endl;
	}
}

// ----------------------------------------------------------------------------

struct SWaveParams
{
	double T;			/** maximal time */
	double dt;			/** time step */
	string fileMesh;
	bool lumping;
};

/** global parameters */
SWaveParams gParams;

// ----------------------------------------------------------------------------

// we're keeping the problem specific functions in a namespace to avoid any ambiguity
namespace wave
{	
	// NOTE: f and g are constant 0 what simplifies the problem

	/** initial value of u */
	double u0(const Vertex& v)
	{
		return exp(-1. * (pow(v.x - 0.1, 2) + pow(v.y - 0.1, 2)) / 0.01);
	}
	
	/** initial value of the time derivative of u */
	double u1(const Vertex& v)
	{
		return 0.;
	}
};

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	clock_t t, tStart, tLoadMesh, tMatA, tMatM, tMatB, tSolve, tEnd;
	double* tSteps;

	tStart = clock();

	// ---
	
	// default values

	gParams.T = 2.25;
	gParams.dt = -1;
	gParams.fileMesh = "data/mesh/cercle1.msh";
	//gParams.fileMesh = "data/mesh/cercle2.msh";
	//gParams.fileMesh = "data/mesh/square_9.msh";
	gParams.lumping = false;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/wave [-lump] -T " << gParams.T << " -dt 0.2 " << gParams.fileMesh << endl;
			cout << "       -lump = use mass lumping for the assembly of M" << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-lump") == 0)
		{
			gParams.lumping = true;
		}
		else if (strcmp(argv[iArg], "-T") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.T;
			}
		}
		else if (strcmp(argv[iArg], "-dt") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.dt;
			}
		}
		else
		{
				gParams.fileMesh = argv[iArg];
		}
	}

	cout << "== Partie II: Wave ==" << endl;
	
	// ----------
	
	// Load the mesh

	RESETCLOCK();
	
	Mesh mesh(gParams.fileMesh.c_str());
	
	CLOCK(tLoadMesh);

	double maxDiameter = mesh.maxIncircleDiameter();
	cout << "Max Diameter: " << maxDiameter << endl;
	
	// ----------
	
	// Prepare time steps
	
	if (gParams.dt <= 0.) 
	{
		gParams.dt = maxDiameter / 3.;
	}
	
	uint nSteps = ceil(gParams.T / gParams.dt);
	tSteps = new double(nSteps);
	
	// ----------
	
	// Assemble the matrices and vectors

	uint dim = mesh.countVertices();

	// @@@
	DUMP(gParams.T);
	DUMP(nSteps);
	DUMP(gParams.dt);
	DUMP(dim);
	
	// vMatV = M
	Sparse uMatU, vMatU, uMatV, M;
	
	{
		// Assemble simple matrices
	
		SparseMap Amap(dim, dim), Mmap(dim, dim), Bmap(dim, dim);
		
		// A
		
		RESETCLOCK();
		
		Amap.constructA(mesh);
		
		CLOCK(tMatA);
		
		// M
		
		RESETCLOCK();
		
		if (gParams.lumping)
		{
			Mmap.constructMlump(mesh);
		}
		else
		{
			Mmap.constructM(mesh);
		}
		
		CLOCK(tMatM);
		
		// B
		
		RESETCLOCK();

		Bmap.constructB(mesh);
		
		CLOCK(tMatB);
	
		// Assemble combined matrices
	
		SparseMap uMatUmap(dim, dim), vMatUmap(dim, dim), uMatVmap(dim, dim);
		
		SparseMap ABmap = Amap;
		ABmap += Bmap;
		
		uMatUmap = ABmap;
		uMatUmap *= -pow(gParams.dt, 2) / 2.;
		uMatUmap += Mmap;
		uMatU = uMatUmap;
		
		vMatUmap = Mmap;
		vMatUmap *= gParams.dt;
		vMatU = vMatUmap;
		
		uMatVmap = ABmap;
		uMatVmap *= -gParams.dt / 2.;
		uMatV = uMatVmap;
		
		// vMatV = M
		M = Mmap;
	}

	// ----------

	// solve the problem
	
	Vector x(dim);
	Vector xLast(dim);
	Vector y(dim);
	Vector yLast(dim);

	Vector originPos(nSteps+1);
	Vector timeVals(nSteps+1);

	
	// initial values
	xLast.constructFunc(mesh, wave::u0);
	yLast.constructFunc(mesh, wave::u1);

	originPos(0) = mesh.eval(0.0, 0.0, xLast);
	timeVals(0) = 0;
	
	for (uint i = 0; i < nSteps; i++)
	{
		cout << "Solving t = " << (i+1) * gParams.dt << "s" << endl;
	
		RESETCLOCK();

		M.newmark(x, y, xLast, yLast, uMatU, vMatU, uMatV);
		
		//DUMP(x.norm2());
		//DUMP(y.norm2());
		
		xLast = x;
		yLast = y;
		
		CLOCK(tSolve);

		// track initial position
		originPos(i+1) = mesh.eval(0.0, 0.0, xLast);
		timeVals(i+1) = (i+1) * gParams.dt;
		
		//tSteps[i] = ((double)tSolve) / CLOCKS_PER_SEC;
		
		stringstream buf;
		buf << "helmholtz_" << (i+1) * gParams.dt;
		string plotFile;
		buf >> plotFile;
		
		PlotMesh plotD(plotFile.c_str(), mesh, x);
		plotD.generate(ePT_GNUPLOT, true, true);
	}
	
	// ----------

	// Evaluate stuff
		
	// stop measuring time, display of graphs is optional
	tEnd = clock();

	// ----------

	// Plot initial data
	
	PlotMesh plotU0("u0", mesh, wave::u0, "Initial distribution");
	plotU0.generate(ePT_GNUPLOT, true);
	
	PlotMesh plotU1("u1", mesh, wave::u1, "Initial velocity");
	plotU1.generate(ePT_GNUPLOT, true);
	
	// Plot last step
	PlotMesh plotX("x", mesh, x, "Final step");
	plotX.generate(ePT_GNUPLOT, true);

	// Plot position of origin over time
	Plot plotOrigin("", timeVals, originPos, "Position of origin over time", "", " w linespoints");
	plotOrigin.generate(ePT_GNUPLOT, true);
	
	/*
	// save the linear system (with our solution for debugging)
	{
		ofstream f("data/linsys/helmholtz.linsys");
		f << AMB << rhs << uh;
	}
	
	// ---
	
	// display computation time

	clock_t tCombinedAssembly = tMatA + tMatM + tAssB + tRhsF + tRhsG;
	clock_t tCombined = tLoadMesh + tSolve + tCombinedAssembly;
	
	// ---
	
	cout << "---" << endl;
	
	LOGPARTTIME("Load Mesh", tLoadMesh, tCombined);
	
	cout << "---" << endl;
	
	LOGPARTTIME("A Assembly", tMatA, tCombined);
	LOGPARTTIME("M Assembly", tMatM, tCombined);
	LOGPARTTIME("B Assembly", tAssB, tCombined);
	
	LOGPARTTIME("F Assembly", tRhsF, tCombined);
	LOGPARTTIME("G Assembly", tRhsG, tCombined);
	
	LOGPARTTIME("Combined Assembly", tCombinedAssembly, tCombined);
	
	cout << "---" << endl;
	
	LOGPARTTIME("Solving", tSolve, tCombined);
	
	cout << "---" << endl;
	
	LOGTIME("Combined time", tCombined);
	
	//cout << "---" << endl;
	//LOGTIME("Total time", tEnd - tStart);
	
	cout << "---" << endl;
	cout << "Computed: " << gParams.fileMesh << " - Nv: " << mesh.countVertices() << 
	", Nt: " << mesh.countTriangles() << ", nE: " << mesh.countEdges() << endl;
	cout << "          Matrix NNZ: " << AMB.sizeNNZ() << " (" << (double)AMB.sizeNNZ() / (AMB.sizeRows() * AMB.sizeColumns()) * 100. << "%)" << endl;
	*/
	
	// ----------
	
	// cleanup
	
	delete[] tSteps;
	
	return 0;
}

// FIN de la partie II :)