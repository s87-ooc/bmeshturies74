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

struct SWaveParams
{
	double T;			/** maximal time */
	uint n;				/** number of partitions of T without 0 */
	string fileMesh;
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
	gParams.n = 10;
	gParams.fileMesh = "data/mesh/cercle1.msh";
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/wave -mesh " << gParams.fileMesh << " -T " << gParams.T << " -n " << gParams.n << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-mesh") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				gParams.fileMesh = argv[iArg];
			}
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
		else if (strcmp(argv[iArg], "-n") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.n;
			}
		}
	}

	cout << "== Partie II: Wave ==" << endl;
	
	// ----------
	
	// Prepare time steps
	double dt = gParams.T / (double)gParams.n;
	tSteps = new double(gParams.n);
	
	// ----------
	
	// Load the mesh

	RESETCLOCK();
	
	Mesh mesh(gParams.fileMesh.c_str());
	
	CLOCK(tLoadMesh);

	// ----------
	
	// Assemble the matrices and vectors

	uint Nv = mesh.countVertices();

	// vMatV = M
	Sparse uMatU, vMatU, uMatV, M;
	
	{
		// Assemble simple matrices
	
		SparseMap Amap(Nv, Nv), Mmap(Nv, Nv), Bmap(Nv, Nv);
		
		// A
		
		RESETCLOCK();
		
		Amap.constructA(mesh);
		
		CLOCK(tMatA);
		
		// M
		
		RESETCLOCK();
		
		Mmap.constructM(mesh);
		
		CLOCK(tMatM);
		
		// B
		
		RESETCLOCK();

		Bmap.constructB(mesh);
		
		CLOCK(tMatB);
	
		// Assemble combined matrices
	
		SparseMap uMatUmap(Nv, Nv), vMatUmap(Nv, Nv), uMatVmap(Nv, Nv);
		
		SparseMap ABmap = Amap;
		ABmap += Bmap;
		
		uMatUmap = ABmap;
		uMatUmap *= -pow(dt, 2) / 2.;
		uMatUmap += Mmap;
		uMatU = uMatUmap;
		
		vMatUmap = Mmap;
		vMatUmap *= dt;
		vMatU = vMatUmap;
		
		uMatVmap = ABmap;
		uMatVmap *= -dt / 2.;
		uMatV = uMatVmap;
		
		// vMatV = M
		M = Mmap;
	}

	// ----------

	// solve the problem
	
	Vector x(Nv);
	Vector xLast(Nv);
	Vector y(Nv);
	Vector yLast(Nv);
	
	// initial values
	xLast.constructFunc(mesh, wave::u0);
	yLast.constructFunc(mesh, wave::u1);
	
	
	for (uint i = 0; i < gParams.n; i++)
	{
		cout << "Solving t = " << (i+1) * dt << "s" << endl;
	
		RESETCLOCK();

		M.newmark(x, y, xLast, yLast, uMatU, vMatU, uMatV);
		
		xLast = x;
		yLast = y;
		
		CLOCK(tSolve);
		
		//tSteps[i] = ((double)tSolve) / CLOCKS_PER_SEC;
		
		// @@@
		PlotMesh plotX("x", mesh, x);
		plotX.generate(true);
	}

	// ----------

	// Evaluate stuff
		
	// stop measuring time, display of graphs is optional
	tEnd = clock();
	
	// ----------

	// Plot initial data
	
	PlotMesh plotU0("u0", mesh, wave::u0);
	plotU0.generate(true);
	
	PlotMesh plotU1("u1", mesh, wave::u1);
	plotU1.generate(true);
	
	// Plot last step
	PlotMesh plotX("x", mesh, x);
	plotX.generate(true);
	
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