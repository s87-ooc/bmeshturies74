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
	uint n;				/** number of partitions of T without 0 */
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
	gParams.n = 43;
	//gParams.n = 4;
	gParams.fileMesh = "data/mesh/cercle1.msh";
	//gParams.fileMesh = "data/mesh/cercle2.msh";
	//gParams.fileMesh = "data/mesh/square_9.msh";
	gParams.lumping = false;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/wave [-lump] -T " << gParams.T << " -n " << gParams.n << " " << gParams.fileMesh << endl;
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
		else
		{
			iArg++;
			
			if (iArg < argc)
			{
				gParams.fileMesh = argv[iArg];
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

	cout << "Max Diameter: " << mesh.maxDiameter() << endl;
	
	// ----------
	
	// Assemble the matrices and vectors

	uint dim = mesh.countVertices();

	// @@@
	DUMP(gParams.T);
	DUMP(gParams.n);
	DUMP(dt);
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
		
		// @@@@
		//_writeMatrix("AB.dat", ABmap);
		
		uMatUmap = ABmap;
		uMatUmap *= -pow(dt, 2) / 2.;
		uMatUmap += Mmap;
		uMatU = uMatUmap;
		
		// @@@@
		//_writeMatrix("uMatU.dat", uMatU);
		
		vMatUmap = Mmap;
		vMatUmap *= dt;
		vMatU = vMatUmap;
		
		// @@@@
		//_writeMatrix("vMatU.dat", vMatU);
		
		uMatVmap = ABmap;
		uMatVmap *= -dt / 2.;
		uMatV = uMatVmap;
		
		// @@@@
		//_writeMatrix("uMatV.dat", uMatV);
		
		// vMatV = M
		M = Mmap;
		
		// @@@@
		//_writeMatrix("M.dat", M);
		
		// @@@
		/*DUMP(-pow(dt, 2) / 2.);
		DUMP(dt);
		DUMP(-dt / 2.);
		Vector vOne(dim);
		for (uint i = 0; i < dim; i++)
		{
			vOne(i) = 1.;
		}
		DUMP((uMatU * vOne).norm2());
		DUMP((vMatU * vOne).norm2());
		DUMP((uMatV * vOne).norm2());
		DUMP((M * vOne).norm2());*/
		
		/*
		{
			ofstream f("uU.txt");
			f << uMatUmap;
		}
		
		{
			ofstream f("vU.txt");
			f << vMatUmap;
		}
		
		{
			ofstream f("uV.txt");
			f << uMatVmap;
		}
		
		{
			ofstream f("vV.txt");
			f << Mmap;
		}
		*/
		return 0;
	}

	// ----------

	// solve the problem
	
	Vector x(dim);
	Vector xLast(dim);
	Vector y(dim);
	Vector yLast(dim);
	
	// initial values
	xLast.constructFunc(mesh, wave::u0);
	yLast.constructFunc(mesh, wave::u1);
	
	// @@@
	DUMP(x.norm2());
	DUMP(xLast.norm2());
	DUMP(y.norm2());
	DUMP(yLast.norm2());
	PlotMesh plot0("x", mesh, xLast, "Initial x");
	plot0.generate(ePT_GNUPLOT, true);
	
	// @@@
	//_writeVector("x0.dat", xLast);
	//_writeVector("y0.dat", yLast);
	
	for (uint i = 0; i < gParams.n; i++)
	{
		cout << "Solving t = " << (i+1) * dt << "s" << endl;
	
		RESETCLOCK();

		M.newmark(x, y, xLast, yLast, uMatU, vMatU, uMatV);
		
		// @@@
		DUMP(x.norm2());
		DUMP(xLast.norm2());
		DUMP(y.norm2());
		DUMP(yLast.norm2());
		
		// @@@
		stringstream buf;
		buf << i+1 << ".dat";
		string fn;
		buf >> fn;
		//_writeVector((string("x")+fn).c_str(), x);
		//_writeVector((string("y")+fn).c_str(), y);
		
		xLast = x;
		yLast = y;
		
		CLOCK(tSolve);
		
		//tSteps[i] = ((double)tSolve) / CLOCKS_PER_SEC;
		
		// @@@
		//PlotMesh plotD("x", mesh, x);
		//plotD.generate(ePT_GNUPLOT, true);
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