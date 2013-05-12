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

//#define DUMPNZ(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { if (fabs(M(i,j)) > 0.0001)	cout << j+1 << " ";	else cout << "-" << " "; } cout << endl; }
//#define DUMPVAL(M) for (uint i = 0; i < M.sizeRows(); i++) {	for (uint j = 0; j < M.sizeColumns(); j++) { cout << M(i,j) << " "; } cout << endl; }

// ----------------------------------------------------------------------------

struct SHelmholtzParams
{
	double kappa;			/** problem parameter */
	double k1;				/** x-component of k */
	double k2;
	string fileMesh;
	bool lumping;
	bool errorBenchmark;
	uint errorCount;
	// ---
	uint* vertexCount;			/** count of vertices for error benchmark */
	uint* triangleCount;		/** count of triangles for error benchmark */
	double* hInner;				/** max inscribed diameter of the triangles */
	double* hOuter;				/** max circumscribed diameter of the triangles */
	double* errors;				/** absolute error on points of the mesh */
	double* errorsL2;			/** L2 error of uh and the projection of the exact solution */
	double* errorsGradientL2;		/** L2 error of the gradient of uh and the projection of the exact solution */
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
	gParams.errorBenchmark = false;
	gParams.errorCount = 1;
	// ---
	gParams.vertexCount = 0;
	gParams.triangleCount = 0;
	gParams.hInner = 0;
	gParams.hOuter = 0;
	gParams.errors = 0;
	gParams.errorsL2 = 0;
	gParams.errorsGradientL2 = 0;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/helmholtz [-lump] [-error n] -kappa " << gParams.kappa << " -k " << gParams.k1 << "," << gParams.k2  << " " << gParams.fileMesh << endl;
			cout << "       -lump = use mass lumping for the assembly of M" << endl;
			cout << "       -error n = do an error benchmark from 0 to n-1 on square and circle meshes" << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-lump") == 0)
		{
			gParams.lumping = true;
		}
		else if (strcmp(argv[iArg], "-error") == 0)
		{
			gParams.errorBenchmark = true;
			
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.errorCount;
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
	
	if (gParams.errorBenchmark)
	{
		gParams.vertexCount = new uint[gParams.errorCount];
		gParams.triangleCount = new uint[gParams.errorCount];
		gParams.hInner = new double[gParams.errorCount];
		gParams.hOuter = new double[gParams.errorCount];
		gParams.errors = new double[gParams.errorCount];
		gParams.errorsL2 = new double[gParams.errorCount];
		gParams.errorsGradientL2 = new double[gParams.errorCount];
	}
	
	// if we're not in benchmark mode, we're just going through once
	
	for (uint iMsh = 0; iMsh < gParams.errorCount; iMsh++)
	{
		
		// selection of the mesh happens automated if in benchmark mode
		if (gParams.errorBenchmark)
		{
			stringstream buf;
			buf << "data/mesh/carre_ffpp_";
			buf << iMsh;
			buf >> gParams.fileMesh;
			
			gParams.fileMesh += ".msh";
		}
		
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
		
		double error = err.norm2();
		double errorL2 = globalL2Error(mesh, u, uh);
		double errorGradL2 = globalL2GradError(mesh, u, uh);
		double hInner = mesh.maxIncircleDiameter();
		double hOuter = mesh.maxCircumcircleDiameter();
		
		if (gParams.errorBenchmark)
		{
			gParams.vertexCount[iMsh] = mesh.countVertices();
			gParams.triangleCount[iMsh] = mesh.countTriangles();
			gParams.hInner[iMsh] = hInner;
			gParams.hOuter[iMsh] = hOuter;
			gParams.errors[iMsh] = error;
			gParams.errorsL2[iMsh] = errorL2;
			gParams.errorsGradientL2[iMsh] = errorGradL2;
		}
		
		cout << "Error: " << error << endl;
		
		cout << "L2 error: " << errorL2 << endl;
		cout << "L2 grad error: " << errorGradL2 << endl;
		cout << "Max inscribed circle diameter: " << hInner << endl;
		cout << "Max circumscribed circle diameter: " << hOuter << endl;

		if (gParams.errorBenchmark)
		{
			cout << "### calculated errors " << iMsh + 1 << " / " << gParams.errorCount << endl << endl;
		}
		else
		{
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
		}

	}

	// ---
	
	// plot error benchmark results
	
	if (gParams.errorBenchmark)
	{
		{
			Vector x(gParams.errorCount);
			Vector y(gParams.errorCount);
			
			for (uint i = 0; i < gParams.errorCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errors[i];
			}
		
			Plot p("helmholtz_hInnerErrors", x, y, "errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
		}
		
		{
			Vector x(gParams.errorCount);
			Vector y(gParams.errorCount);
			
			for (uint i = 0; i < gParams.errorCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsL2[i];
			}
		
			Plot p("helmholtz_hInnerErrorsL2", x, y, "L2 errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
		}
		
		{
			Vector x(gParams.errorCount);
			Vector y(gParams.errorCount);
			
			for (uint i = 0; i < gParams.errorCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsGradientL2[i];
			}
		
			Plot p("helmholtz_hInnerErrorsGradientL2", x, y, "L2 grad errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
		}
	}
	
	// ---
	
	// cleanup
	
	SAFE_ARRDELETE(gParams.vertexCount);
	SAFE_ARRDELETE(gParams.triangleCount);
	SAFE_ARRDELETE(gParams.hInner);
	SAFE_ARRDELETE(gParams.hOuter);
	SAFE_ARRDELETE(gParams.errors);
	SAFE_ARRDELETE(gParams.errorsL2);
	SAFE_ARRDELETE(gParams.errorsGradientL2);
	
	return 0;
}

// FIN de la partie I :)