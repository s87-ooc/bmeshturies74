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

using namespace std;

// ----------------------------------------------------------------------------

#define RESETCLOCK _t = clock
#define CLOCK(T) T = clock() - _t
#define LOGTIME(S, T)cout << S << ": " << ((float)T) / CLOCKS_PER_SEC << "s (" << T << " clocks)" << endl;
#define LOGPARTTIME(S, T, TMAX)cout << S << ": " << ((float)T) / CLOCKS_PER_SEC << "s (" << T << " clocks) " << setprecision(2) << fixed << (double)T / TMAX * 100. << "%" << endl;

clock_t _t;

#define DUMPNZ(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { if (fabs(M(i,j)) > 0.0001)	cout << j+1 << " ";	else cout << "-" << " "; } cout << endl; }

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

namespace helmholtz
{

	double f(const Vertex& v)
	{
		return 0.;
	}

	double g(const Vertex& v)
	{
		// make sure we have a vertex on the edge
		if (v.label != 1)
		{
			return 0.;
		}

		// constant is the same for any border
		double kapxk = gParams.kappa * (gParams.k1 * v.x + gParams.k2 * v.y);
		
		double n1, n2;
		
		// find the proper edge if we have a square
		double tol = 0.00000001;
		if (fabs(v.x) < tol)
		{
			if (fabs(v.y) < tol || fabs(v.y - 1.) < tol)
				return 0.;
			
			// left border
			n1 = -1.;
			n2 = 0.;
		}
		else if (fabs(v.x - 1.) < tol)
		{
			if (fabs(v.y) < tol || fabs(v.y - 1.) < tol)
				return 0.;
		
			// right border
			n1 = 1.;
			n2 = 0.;
		}
		else if (fabs(v.y) < tol)
		{
			// bottom border
			n1 = 0.;
			n2 = -1.;
		}
		else if (fabs(v.y - 1.) < tol)
		{
			// top border
			n1 = 0.;
			n2 = 1.;
		}

		double kn = gParams.k1 * n1 + gParams.k2 * n2;
		
		double val = sin(kapxk) + gParams.kappa * cos(kapxk) * kn;
		//cout << v.x << " " << v.y << " " << val << endl;
		
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
	
	// get filenames and parameters from the console
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: helmholtz -mesh " << gParams.fileMesh << " -kappa " << gParams.kappa << " -k " << gParams.k1 << " " << gParams.k2 << endl;
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
	
	Mesh m(gParams.fileMesh.c_str());
	
	CLOCK(tLoadMesh);

	// ----------
	
	// Assembe the matrices
	
	uint nV = m.countVertices();
	SparseMap Amap(nV, nV), Mmap(nV, nV), Bmap(nV, nV);
	
	RESETCLOCK();
	
	// @@@
	Amap.constructA(m);
	
	CLOCK(tMatA);

	if (0)
	{
		// Verification
		cout << endl << "A" << endl << endl;
		Vector k(Amap.sizeRows());
		for (uint i = 0; i < k.size(); i++)
		{
			k(i) = 1.;
		}
		Sparse _A(Amap);
		Vector Ak = _A*k;
		cout << "Ak: " << Ak.norm2() << endl;
		
		const SparseMap& refA = Amap;
		DUMPNZ(refA);
		
		cout << endl;
		
		DUMPNZ(_A);
	
	}
	
	// ---
	
	RESETCLOCK();
	
	Mmap.constructM(m);
	
	CLOCK(tMatM);
	
	if (0)
	{	
		cout << endl << "M" << endl << endl;
	
		const SparseMap& refM = Mmap;
		DUMPNZ(refM);
	
		cout << endl;
		double area = 0.;
		for (uint t = 0; t < m.countTriangles(); t++)
		{
			area += m.T[t].area;
		}
		cout << "Area: " << area << endl;
		
		cout << "Sum over Mij: ";
		double areaM = 0.;
		Sparse M(Mmap);
		for (uint i = 0; i < M.sizeRows(); i++)
		{
			for (uint j = 0; j < M.sizeColumns(); j++)
			{
				areaM += M(i,j);
			}
		}
		cout << areaM << endl;
		
		DUMPNZ(M);
	}
	
	// ---
	
	RESETCLOCK();
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// @@@
	Bmap.constructB(m);
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	CLOCK(tAssB);

	cout << endl << "B" << endl << endl;
	
	const SparseMap& refB = Bmap;
	//DUMPNZ(refB);
	
	cout << endl;
	
	Sparse B(Bmap);
	//DUMPNZ(B);
	
	// ----------	
	
	// Assemble rhs
	
	Vector rhsF(nV), rhsG(nV);
	
	RESETCLOCK();
	
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	rhsF.constructFunc(m, &helmholtz::f);
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	CLOCK(tRhsF);
	
	// ---
	
	RESETCLOCK();
	
	// @@@
	rhsG.constructFuncIntSurf(m, &helmholtz::g);
	
	CLOCK(tRhsG);
	
	// ----------
	
	// solve the stationary problem

	cout << "A1" << endl;
	//DUMPNZ(refA);
	
	// @@@
	Mmap *= -pow(gParams.kappa, 2);	
	
	cout << "M2" << endl;
	//DUMPNZ(Mmap);
	
	// @@@
	Amap += Mmap;
	
	cout << "A2" << endl;
	//DUMPNZ(refA);
	
	// @@@
	Amap += Bmap;
	
	cout << "A3" << endl;
	//DUMPNZ(refA);
	
	Sparse AMB(Amap);
	cout << "ASparse" << endl;
	//DUMPNZ(AMB);
	
	//Vector rhs = rhsF + rhsG;
	
	//return 0;
	
	RESETCLOCK();

	Vector solution = AMB.conjGradient(rhsG);
	//Vector solution = AMB.jacobi(rhsG);
	
	CLOCK(tSolve);
	
	// ----------
	
	// Evaluate stuff
	
	// stop measuring time, display of graphs is optional
	tEnd = clock();
	
	// ----------
	
	// plot rhs functions
	
	//PlotMesh plotF("f", &m, &helmholtz::f);
	//plotF.generate(true);
	
	PlotMesh plotG("g", &m, &helmholtz::g);
	//plotG.generate(true);

	// plot our solution

	PlotMesh plotUh("uh", &m, &solution);
	plotUh.generate(true);
	
	// plot exact solution
	
	PlotMesh plotU("u", &m, &helmholtz::u);
	plotU.generate(true);
	
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
	cout << "Computed: " << gParams.fileMesh << " - Nv: " << m.countVertices() << 
	", Nt: " << m.countTriangles() << ", nE: " << m.countEdges() << endl;
	
	return 0;
}