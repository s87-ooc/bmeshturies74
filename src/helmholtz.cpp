/*************************************************************

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

#define DUMPNZ(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { if (fabs(M(i,j)) > 0.0001)	cout << j+1 << " ";	else cout << "-" << " "; } cout << endl; }
#define DUMPVAL(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { cout << M(i,j) << " "; } cout << endl; }

// ----------------------------------------------------------------------------

struct SHelmholtzParams
{
	double kappa;
	double k1;
	double k2;
	string fileMesh;
};

/** global parameters */
SHelmholtzParams gParams;

// ----------------------------------------------------------------------------

// we're keeping the problem specific functions in a namespace to avoid any ambiguity
namespace helmholtz
{
	// f is just
	double f(const Vertex& v)
	{
		return 0.;
	}

	// TODO: what is the normal on the corners of the square?
	double g(const Vertex& v)
	{
		// make sure we have a vertex on the edge
		if (v.label != 1)
		{
			return 0.;
		}

		// this constant is the same for any border
		double kapxk = gParams.kappa * (gParams.k1 * v.x + gParams.k2 * v.y);
		
		// normal vector n = (n1, n2)
		double n1, n2;
		
		// TODO: toggle corner point treatment here to see the impact on the solution
		bool checkCorners = true;
		
		// determine the part of the border the point belongs to
		double tol = 0.00000001;

		// LEFT BORDER: n = (-1, 0) on the edge
		if (fabs(v.x) < tol)
		{
			if (checkCorners)
			{
				// LOWER LEFT CORNER
				if (fabs(v.y) < tol)
				{
					n1 = -1. / sqrt(2.);
					n2 = -1. / sqrt(2.);
				}
				// UPPER LEFT CORNER
				else if (fabs(v.y - 1.) < tol)
				{
					n1 = -1. / sqrt(2.);
					n2 = 1. / sqrt(2.);
				}
				else
				{
					n1 = -1.;
					n2 = 0.;
				}
			}
			else
			{
				n1 = -1.;
				n2 = 0.;
			}
		}
		// RIGHT BORDER: n = (1, 0) on the edge
		else if (fabs(v.x - 1.) < tol)
		{
			if (checkCorners)
			{
				// LOWER RIGHT CORNER
				if (fabs(v.y) < tol)
				{
					n1 = -1 / sqrt(2.);
					n2 = -1 / sqrt(2.);
				}
				// UPPER RIGHT CORNER
				else if (fabs(v.y - 1.) < tol)
				{
					n1 = -1 / sqrt(2.);
					n2 = -1 / sqrt(2.);
				}
				else
				{
					n1 = 1.;
					n2 = 0.;
				}
			}
			else
			{
				// right border
				n1 = 1.;
				n2 = 0.;
			}
		}
		// if we had a corner it was already treated above
		// BOTTOM BORDER: n = (0, -1) on the edge
		else if (fabs(v.y) < tol)
		{
			n1 = 0.;
			n2 = -1.;
		}
		// TOP BORDER: n = (0, 1) on the edge
		else if (fabs(v.y - 1.) < tol)
		{
			n1 = 0.;
			n2 = 1.;
		}

		double kn = gParams.k1 * n1 + gParams.k2 * n2;
		
		double val = sin(kapxk) + gParams.kappa * cos(kapxk) * kn;
		//cout << v.x << " " << v.y << " " << val << endl;
		
		return val;
	}

	double g2(const Vertex& v, const BoundEdge& e)
	{
		double kapxk = gParams.kappa * (gParams.k1 * v.x + gParams.k2 * v.y);

		Vector n(2);
		n.constructNormal(e);

		double kn = n(0)*gParams.k1 + n(1)*gParams.k2;
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
	clock_t t, tStart, tLoadMesh, tMatA, tMatM, tAssB, tRhsF, tRhsG, tSolve, tEnd;
	
	tStart = clock();

	// ---
	
	// default values

	gParams.kappa = 10.;
	gParams.k1 = 1./sqrt(2.);
	gParams.k2 = 1./sqrt(2.);
	gParams.fileMesh = "data/mesh/carre1.msh";
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/helmholtz -mesh " << gParams.fileMesh << " -kappa " << gParams.kappa << " -k " << gParams.k1 << " " << gParams.k2 << endl;
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
		// TODO: k as parameter
	}

	cout << "== Partie I: Helmholtz ==" << endl;

	// ----------
	
	// Load the mesh

	RESETCLOCK();
	
	Mesh mesh(gParams.fileMesh.c_str());
	
	CLOCK(tLoadMesh);

	// ----------
	
	// Assemble the matrices and vectors
	
	uint Nv = mesh.countVertices();
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
	
	CLOCK(tAssB);
	
	// F (rhs)
	
	Vector rhsF(Nv), rhsG(Nv);
	
	RESETCLOCK();
	
	rhsF.constructFuncIntSurf(mesh, helmholtz::f);
	
	CLOCK(tRhsF);
	
	// G (rhs)
	
	RESETCLOCK();
	
	//rhsG.constructFuncIntSurf(mesh, helmholtz::g);
	rhsG.constructFuncSurf(mesh, helmholtz::g2);
	
	CLOCK(tRhsG);
	
	// ----------
	
	// verification of constructed matrices according to the course
	
	// dumping may make sense for "bin/helmholtz -mesh data/mesh/square_9.msh"
	
	bool verify = false;
	bool dump = false;
	
	// A

	if (verify || dump)
		cout << endl << "A" << endl << endl;
		
	if (dump)
	{
		DUMPNZ(Amap);
		cout << endl;
		DUMPVAL(Amap);
		cout << endl;
	}
	
	if (verify)
	{
		// TEST: A * (1, ..., 1) = 0
		
		Vector k(Amap.sizeRows());
		for (uint i = 0; i < k.size(); i++)
		{
			k(i) = 1.;
		}
		Sparse _A(Amap);
		Vector Ak = _A*k;
		cout << "Ak: " << Ak.norm2() << endl;
	}
	
	// M

	if (verify || dump)
		cout << endl << "M" << endl << endl;
		
	if (dump)
	{
		DUMPNZ(Mmap);
		cout << endl;
		DUMPVAL(Mmap);
		cout << endl;
	}
	
	if (verify)
	{
		// TEST: Sum Mij = area of mesh		

		double area = 0.;
		for (uint t = 0; t < mesh.countTriangles(); t++)
		{
			area += mesh.T[t].area;
		}
		cout << "Area: " << area << endl;
		
		cout << "Sum over Mij: ";
		double areaM = 0.;
		for (uint i = 0; i < Mmap.sizeRows(); i++)
		{
			for (uint j = 0; j < Mmap.sizeColumns(); j++)
			{
				areaM += Mmap(i,j);
			}
		}
		cout << areaM << endl;
	}
	
	// B
	
	if (verify || dump)
		cout << endl << "B" << endl << endl;
		
	if (dump)
	{
		DUMPNZ(Bmap);
		cout << endl;
		DUMPVAL(Bmap);
		cout << endl;
	}
	
	if (verify)
	{
	}
	
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
	
	CLOCK(tSolve);
	
	// ----------
	
	// Evaluate stuff
	
	// calculate exact solution
	Vector u(Nv);
	u.constructFunc(mesh, helmholtz::u);
	
	// calculate error
	Vector err = u - uh;
	
	// stop measuring time, display of graphs is optional
	tEnd = clock();
	
	cout << "Error: " << err.norm2() << endl;
	
	// ----------
	
	// Plot problem functions and results
	
	PlotMesh plotF("f", mesh, helmholtz::f);
	plotF.generate(true);
	
	// TODO: check the impact of g to the solution
	PlotMesh plotG("g", mesh, helmholtz::g);
	plotG.generate(true);

	// our solution

	PlotMesh plotUh("uh", mesh, uh);
	plotUh.generate(true);
	
	// exact solution
	
	PlotMesh plotU("u", mesh, u);
	plotU.generate(true);
	
	// error
	
	PlotMesh plotErr("err", mesh, err);
	plotErr.generate(true);
	
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
	
	return 0;
}

// FIN de la partie I :)