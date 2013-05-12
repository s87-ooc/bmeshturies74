﻿/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: helmholtz.cpp

 Description: solution de la première partie

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

// dump info on non-zero entries in (small!) matrices

//#define DUMPNZ(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { if (fabs(M(i,j)) > 0.0001)	cout << j+1 << " ";	else cout << "-" << " "; } cout << endl; }
//#define DUMPVAL(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { cout << M(i,j) << " "; } cout << endl; }

// ----------------------------------------------------------------------------

struct SHelmholtzParams
{
	double kappa;
	double k1;
	double k2;
	string fileMesh;
	bool lumping;
};

/** global parameters */
SHelmholtzParams gParams;

// ----------------------------------------------------------------------------

// we're keeping the problem specific functions in a namespace to avoid any ambiguity
namespace helmholtz
{
	// f is just 0
	double f(const Vertex& v)
	{
		return 0.;
	}

	double g(const Vertex& v, const BoundEdge& e)
	{
		double kapxk = gParams.kappa * (gParams.k1 * v.x + gParams.k2 * v.y);

		Vector n(2);
		n.constructNormal(e);

		double kn = n(0) * gParams.k1 + n(1) * gParams.k2;
		double val = sin(kapxk) + gParams.kappa * cos(kapxk) * kn;

		return val;
	}
	
	/** exact solution */
	double u(const Vertex& v)
	{
		return sin(gParams.kappa * (gParams.k1 * v.x + gParams.k2 * v.y));
	}

};

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	clock_t t, tStart, tLoadMesh, tMatA, tMatM, tMatB, tRhsF, tRhsG, tSolve, tEnd;
	
	tStart = clock();

	// ---
	
	// default values

	gParams.kappa = 10.;
	gParams.k1 = 1./sqrt(2.);
	gParams.k2 = 1./sqrt(2.);
	gParams.fileMesh = "data/mesh/carre1.msh";
	gParams.lumping = false;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/helmholtz [-lump] -kappa " << gParams.kappa << " -k " << gParams.k1 << "," << gParams.k2  << " " << gParams.fileMesh << endl;
			cout << "       -lump = use mass lumping for the assembly of M" << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-lump") == 0)
		{
			gParams.lumping = true;
		}
		else if (strcmp(argv[iArg], "-kappa") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.kappa;
			}
		}
		else if (strcmp(argv[iArg], "-k") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				char* sDim = strtok(argv[iArg], ",");
				
				{
					stringstream buf;
					buf << sDim;
					buf >> gParams.k1;
				}
				
				sDim = strtok(0, ",");
				
				{
					stringstream buf;
					buf << sDim;
					buf >> gParams.k2;
				}
			}
			
			cout << "k1 " << gParams.k1 << " k2 " << gParams.k2 << endl;
		}
		else
		{
			gParams.fileMesh = argv[iArg];
		}
	}

	cout << "== Partie I: Helmholtz ==" << endl;

	// ----------

	{
		ifstream fileTest(gParams.fileMesh.c_str());
		
		if (!fileTest)
		{
			cout << "No valid mesh " << gParams.fileMesh << endl;
			return 0;
		}
	}
	
	// Load the mesh

	RESETCLOCK();
	
	Mesh mesh(gParams.fileMesh.c_str());
	
	CLOCK(tLoadMesh);
	
	// ----------
	
	// Assemble the matrices and vectors
	
	uint dim = mesh.countVertices();
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
	
	// F (rhs)
	
	Vector rhsF(dim), rhsG(dim);
	
	RESETCLOCK();
	
	rhsF.constructFuncIntSurf(mesh, helmholtz::f);
	
	CLOCK(tRhsF);
	
	// G (rhs)
	
	RESETCLOCK();

	rhsG.constructFuncSurf(mesh, helmholtz::g);
	
	CLOCK(tRhsG);
	
	// ----------

	// construct the linear system

	Mmap *= -pow(gParams.kappa, 2);	
	Amap += Mmap;
	Amap += Bmap;

	// AMB = A - kappa^2 * M + B
	// convert the matrix to Sparse so we can apply a solver
	Sparse AMB(Amap);

	Vector rhs = rhsF + rhsG;
	
	// ----------

	// solve the stationary problem	
	
	RESETCLOCK();

	Vector uh = AMB.conjGradient(rhs);
	//Vector uh = AMB.jacobi(rhs);
	//Vector uh = LU.jacobi(rhs);
	
	CLOCK(tSolve);
	
	// ----------
	
	// Evaluate stuff
	
	// calculate exact solution
	Vector u(dim);
	u.constructFunc(mesh, helmholtz::u);
	
	// calculate error
	Vector err = u - uh;
	
	// stop measuring time, display of graphs is optional
	tEnd = clock();
	
	cout << "Error: " << err.norm2() << endl;
	
	cout << "L2 error: " << globalL2Error(mesh, u, uh) << endl;
	cout << "L2 grad error: " << globalL2GradError(mesh, u, uh) << endl;
	cout << "Max circumscribed circle diameter: " << mesh.maxCircumcircleDiameter() << endl;
	cout << "Max inscribed circle diameter: " << mesh.maxIncircleDiameter() << endl;
	// ----------
	
	// boundary conditions
	
	//PlotMesh plotF("f", mesh, helmholtz::f);
	//plotF.generate(ePT_GNUPLOT, true);
	
	//PlotMesh plotG("g", mesh, helmholtz::g);
	//plotG.generate(ePT_GNUPLOT, true);
	
	//PlotMesh plotG("g", mesh, helmholtz::g);
	//plotG.generate(ePT_GNUPLOT, true);

	// our solution

	PlotMesh plotUh("helmholtz_uh", mesh, uh, "Solution FEM");
	plotUh.generate(ePT_MEDIT);
	plotUh.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
	
	// exact solution
	
	PlotMesh plotU("helmholtz_u", mesh, u, "Solution Exacte");
	plotU.generate(ePT_MEDIT);
	plotU.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
	
	// error
	
	PlotMesh plotErr("helmholtz_err", mesh, err, "Erreur");
	plotErr.generate(ePT_MEDIT);
	plotErr.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
	
	// ---
	
	// display computation time

	clock_t tCombinedAssembly = tMatA + tMatM + tMatB + tRhsF + tRhsG;
	clock_t tCombined = tLoadMesh + tSolve + tCombinedAssembly;
	
	// ---
	
	cout << "---" << endl;
	
	LOGPARTTIME("Load Mesh", tLoadMesh, tCombined);
	
	cout << "---" << endl;
	
	LOGPARTTIME("A Assembly", tMatA, tCombined);
	LOGPARTTIME("M Assembly", tMatM, tCombined);
	LOGPARTTIME("B Assembly", tMatB, tCombined);
	
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
	
	return 0;
}

// FIN de la partie I :)