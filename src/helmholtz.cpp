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

struct SHelmholtzParams
{
	double kappa;				/** problem parameter */
	double k1;					/** x-component of k */
	double k2;					/** y-component of k */
	string fileMesh;
	bool lumping;
	uint meshCount;
	bool save;					/** whether to save the plot */
	bool test;					/** performs tests on precision and computation time */
	// values below are used in test mode, [0] = normal [1] = mass lumping:
	uint* vertexCount;				/** count of vertices for test mode */
	uint* triangleCount;			/** count of triangles for test mode */
	double* hInner;					/** max inscribed diameter of the triangles */
	double* hOuter;					/** max circumscribed diameter of the triangles */
	double* errors[2];				/** absolute error on points of the mesh */
	double* errorsL2[2];			/** L2 error of uh and the projection of the exact solution */
	double* errorsGradientL2[2];	/** L2 error of the gradient of uh and the projection of the exact solution */
	double* solveTime[2];			/** time needed for solving the system */
	double* solveAssTime[2];		/** time needed for assembling and solving the system */
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

	// commercial break
	
	cout << "  -- Helmholtz (Partie 1), Projet Final [MM031], 2013" << endl;
	cout << "     Copyright (C) K. Podkanski, S. Stamenkovic" << endl;
	
	// ---
	
	// default values

	gParams.kappa = 10.;
	gParams.k1 = 1./sqrt(2.);
	gParams.k2 = 1./sqrt(2.);
	gParams.fileMesh = "";
	gParams.lumping = false;
	gParams.meshCount = 4;
	gParams.save = false;
	gParams.test = false;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/helmholtz [-lump] [-save] [-test] -kappa " << gParams.kappa << " -k " << gParams.k1 << "," << gParams.k2  << " " << gParams.fileMesh << endl;
			cout << "       -lump = use mass lumping for the assembly of M" << endl;
			cout << "       -save = generate PNGs of the plots" << endl;
			cout << "       -test = perform some tests on different square meshes" << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-lump") == 0)
		{
			gParams.lumping = true;
		}
		else if (strcmp(argv[iArg], "-save") == 0)
		{
			gParams.save = true;
		}
		else if (strcmp(argv[iArg], "-test") == 0)
		{
			gParams.test = true;
			
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.meshCount;
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
			
			Vector k(2);
			k(0) = gParams.k1;
			k(1) = gParams.k2;
			
			// normalize
			gParams.k1 *= 1. / k.norm2();
			gParams.k2 *= 1. / k.norm2();
		}
		else
		{
			gParams.fileMesh = argv[iArg];
			gParams.meshCount = 1;
		}
	}

	// ----------
	
	gParams.vertexCount = new uint[gParams.meshCount];
	gParams.triangleCount = new uint[gParams.meshCount];
	gParams.hInner = new double[gParams.meshCount];
	gParams.hOuter = new double[gParams.meshCount];
	gParams.errors[0] = new double[gParams.meshCount];
	gParams.errors[1] = new double[gParams.meshCount];
	gParams.errorsL2[0] = new double[gParams.meshCount];
	gParams.errorsL2[1] = new double[gParams.meshCount];
	gParams.errorsGradientL2[0] = new double[gParams.meshCount];
	gParams.errorsGradientL2[1] = new double[gParams.meshCount];
	gParams.solveTime[0] = new double[gParams.meshCount];
	gParams.solveTime[1] = new double[gParams.meshCount];
	gParams.solveAssTime[0] = new double[gParams.meshCount];
	gParams.solveAssTime[1] = new double[gParams.meshCount];

	// ---
	
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	
	// iterate through meshes we have to
	
	for (int iMsh = 0; iMsh < gParams.meshCount; iMsh++)
	{
	
		// selection of the mesh happens automated if in test mode
		if (gParams.test)
		{
			stringstream buf;
			buf << "data/mesh/carre_ffpp_";
			buf << iMsh;
			buf >> gParams.fileMesh;
			
			gParams.fileMesh += ".msh";
		}
		else if (gParams.meshCount > 1)
		{
			stringstream buf;
			buf << "data/mesh/carre";
			buf << iMsh + 1;
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

		cout << endl << "Loading mesh..." << endl;
		
		// Load the mesh

		RESETCLOCK();
		
		Mesh mesh(gParams.fileMesh.c_str());
		
		CLOCK(tLoadMesh);
		
		// ----------

		cout << "Solving the linear system..." << endl;
		
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

		// convert back to sparselil for use with lu
		//SparseLIL AMBlil(AMB,1);


		Vector rhs = rhsF + rhsG;
		
		// ----------

		// solve the stationary problem	
		
		RESETCLOCK();

		Vector uh = AMB.conjGradient(rhs);
		//Vector uh = AMB.jacobi(rhs);
		//Vector uh = AMBlil.LU(rhs);
		
		CLOCK(tSolve);
		
		// ----------
		
		// Evaluate stuff
		
		// calculate exact solution
		Vector u(dim);
		u.constructFunc(mesh, helmholtz::u);
		
		// calculate error
		Vector err = u - uh;

		// ---
		
		// combine clock times
		
		// display computation time

		clock_t tCombinedAssembly = tMatA + tMatM + tMatB + tRhsF + tRhsG;
		clock_t tCombined = tSolve + tCombinedAssembly;
		
		// calculate errors
		
		double error = err.norm2();
		double errorL2 = globalL2Error(mesh, u, uh);
		double errorGradL2 = globalL2GradError(mesh, u, uh);
		double hInner = mesh.maxIncircleDiameter();
		double hOuter = mesh.maxCircumcircleDiameter();

		gParams.vertexCount[iMsh] = mesh.countVertices();
		gParams.triangleCount[iMsh] = mesh.countTriangles();
		gParams.hInner[iMsh] = hInner;
		gParams.hOuter[iMsh] = hOuter;
		int idx = gParams.lumping ? 1 : 0;
		gParams.errors[idx][iMsh] = error;
		gParams.errorsL2[idx][iMsh] = errorL2;
		gParams.errorsGradientL2[idx][iMsh] = errorGradL2;
		gParams.solveTime[idx][iMsh] = ((double)tSolve)/CLOCKS_PER_SEC;
		gParams.solveAssTime[idx][iMsh] = ((double)tCombined)/CLOCKS_PER_SEC;

		cout << "---" << endl;
		
		LOGTIME("Load Mesh", tLoadMesh);
		
		cout << "---" << endl;
		
		LOGPARTTIME("A Assembly", tMatA, tCombinedAssembly);
		LOGPARTTIME("M Assembly", tMatM, tCombinedAssembly);
		LOGPARTTIME("B Assembly", tMatB, tCombinedAssembly);
		
		//LOGPARTTIME("F Assembly", tRhsF, tCombinedAssembly);
		//LOGPARTTIME("G Assembly", tRhsG, tCombinedAssembly);
		
		LOGPARTTIME("Combined Assembly", tCombinedAssembly, tCombined);
		
		cout << "---" << endl;
		
		LOGTIME("Solving", tSolve);
		
		cout << "---" << endl;
		
		LOGTIME("Combined time", tCombined);
		
		cout << "---" << endl;
		cout << "Computed: " << gParams.fileMesh << " - Nv: " << mesh.countVertices() << 
		", Nt: " << mesh.countTriangles() << ", nE: " << mesh.countEdges() << ", h: " << hInner << endl;
		cout << "          Matrix NNZ: " << AMB.sizeNNZ() << " (" << (double)AMB.sizeNNZ() / (AMB.sizeRows() * AMB.sizeColumns()) * 100. << "%)" << endl;

		cout << "---" << endl;
		
		cout << "Error: " << error << endl;
		
		cout << "L2 error: " << errorL2 << endl;
		cout << "L2 grad error: " << errorGradL2 << endl;
		//cout << "Max inscribed circle diameter: " << hInner << endl;
		//cout << "Max circumscribed circle diameter: " << hOuter << endl;
		
		// ---

		// plot initial data only if we're having only one mesh to check
		
		if (gParams.meshCount == 1)
		{
			// boundary conditions

			//PlotMesh plotF("helmholtz_f", mesh, helmholtz::f);
			//plotF.generate(ePT_GNUPLOT, true);

			//PlotMesh plotG("helmholtz_g", mesh, helmholtz::g);
			//plotG.generate(ePT_GNUPLOT, true);

			// exact solution

			PlotMesh plotU("p1t2_u", mesh, u, "Solution Exacte");
			plotU.generate(ePT_MEDIT);
			plotU.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
			if (gParams.save) { plotU.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }
		}
		
		if(!gParams.test && !gParams.lumping)
		{
			ostringstream buf;
			buf	<< iMsh + 1;
			
			// our solution

			PlotMesh plotUh(("p1t2_uh_" + buf.str()).c_str(), mesh, uh, "Solution FEM");
			plotUh.generate(ePT_MEDIT);
			plotUh.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
			if (gParams.save) { plotUh.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }
			
			// error
			
			PlotMesh plotErr(("p1t2_err_" + buf.str()).c_str(), mesh, err, "Erreur");
			plotErr.generate(ePT_MEDIT);
			plotErr.generate(ePT_GNUPLOT_SURF, true, false, "", "", 20);
			if (gParams.save) { plotErr.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }
		}		
		
		// ---

		if (iMsh == gParams.meshCount - 1 && !gParams.lumping)
		{
			// reset the iteration and calculate values with mass lumping
			cout << endl << "### Calculate values with mass lumping" << endl << endl;
			iMsh = -1;
			gParams.lumping = true;
		}
		
	}

	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	
	// plot combined data
	
	if (gParams.meshCount > 1)
	{
		/*{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errors[0][i];
			}
		
			Plot p("helmholtz_hInnerErrors", x, y, "errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}*/
		
		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsL2[0][i];
			}
		
			string pName = "p1t3_errL2";
			pName += gParams.test ? "GEN" : "";
			Plot p(pName.c_str(), x, y, "L2 errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}
		
		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsGradientL2[0][i];
			}
			
			string pName = "p1t3_errgradL2";
			pName += gParams.test ? "GEN" : "";
			Plot p(pName.c_str(), x, y, "L2 grad errors over inscribed circle", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}
		
		// TODO
		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			Vector ylump(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsL2[0][i];
			}
		
			Plot p("p2t4_precision", x, y, "Mass Lumping: Precision", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}
		
		// TODO
		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			Vector ylump(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.errorsGradientL2[0][i];
			}
		
			Plot p("p2t4_precisionGrad", x, y, "Mass Lumping: Precision Gradient", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}

		// TODO
		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			Vector ylump(gParams.meshCount);
			
			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.solveTime[0][i];
			}
		
			Plot p("p2t4_time1", x, y, "Mass Lumping: Time", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}
		
	// ---
	
	// cleanup
	
	SAFE_ARRDELETE(gParams.vertexCount);
	SAFE_ARRDELETE(gParams.triangleCount);
	SAFE_ARRDELETE(gParams.hInner);
	SAFE_ARRDELETE(gParams.hOuter);
	SAFE_ARRDELETE(gParams.errors[0]);
	SAFE_ARRDELETE(gParams.errors[1]);
	SAFE_ARRDELETE(gParams.errorsL2[0]);
	SAFE_ARRDELETE(gParams.errorsL2[1]);
	SAFE_ARRDELETE(gParams.errorsGradientL2[0]);
	SAFE_ARRDELETE(gParams.errorsGradientL2[1]);
	SAFE_ARRDELETE(gParams.solveTime[0]);
	SAFE_ARRDELETE(gParams.solveTime[1]);
	SAFE_ARRDELETE(gParams.solveAssTime[0]);
	SAFE_ARRDELETE(gParams.solveAssTime[1]);

	tEnd = clock();
	
	cout << "---" << endl;
	LOGTIME("Total time", tEnd - tStart);
	
	return 0;
}

// FIN de la partie I :)