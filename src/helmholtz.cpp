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

enum EHelmholtzMode
{
	eHM_SIMPLE = 0,
	eHM_TIMETEST,
	eHM_PRECISIONTEST,
	eHM_COMBINEDTEST
};

enum ETime
{
	eT_START = 0,
	eT_LOADMESH,
	eT_MATA,
	eT_MATM,
	eT_MATB,
	eT_RHSF,
	eT_RHSG,
	eT_SOLVE,
	eT_END
};

enum ESolver
{
	eS_CG = 0,
	eS_JACOBI,
	eS_LU
};

// ----------------------------------------------------------------------------

/** parse the commandline parameters */
int parseCmd(int argc, char* argv[]);

// ----------

/** simple run, calculate the solution over the specified mesh(es) and visualize the error with the exact solution */
int simple();

/** timetest, solve the problem on the meshes and compare calculation times, optionally with lumping */
int timetest();

/** precisiontest, solve the problem over the meshes and compare precision, optionally with lumping */
int precisiontest();

// ----------

/** print stats after solving a problem */
void printStats(const char* meshFile, const Mesh& mesh, clock_t times[],
				double error, double errorL2, double errorGradientL2);

// ----------------------------------------------------------------------------

struct SHelmholtzParams
{
	EHelmholtzMode mode;
	std::vector<string> files;	/** list of mesh files we should perform our calculations on */
	string outPath;				/** path where generated plot files, etc should be stored in */
	bool calcMDefault;			/** whether to assemble A using the default method, can be combined with calcMLumping */
	bool calcMLumping;			/** whether to assemble A using mass lumping, can be combined with calcMDefault */
	ESolver solver;				/** which solver should be used for the linear system(s) */
	bool quiet;					/** if true, gnuplot won't be launched to display plots */
	bool generatePNG;			/** if true, gnuplot will be launched to generate a png */
	// ---
	double kappa;				/** problem parameter */
	double k1;					/** x-component of k */
	double k2;					/** y-component of k */
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
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	// commercial break
	
	cout << "  -- Helmholtz (Partie 1), Projet Final [MM031], 2013" << endl;
	cout << "     Copyright (C) K. Podkanski, S. Stamenkovic" << endl << endl;

	// ---
	
	// default values

	gParams.mode = eHM_SIMPLE;
	gParams.files.push_back("data/mesh/carre1.msh");
	gParams.outPath = "data/plots/";
	gParams.calcMDefault = false;	// will be set to true after parsing if no other method specified
	gParams.calcMLumping = false;
	gParams.solver = eS_CG;
	gParams.quiet = false;
	gParams.generatePNG = false;
	// ---
	gParams.kappa = 10.;
	gParams.k1 = 1./sqrt(2.);
	gParams.k2 = 1./sqrt(2.);

	// ---

	// parse cmdline arguments, gives -1 if we should continue the program

	int iReturn = parseCmd(argc, argv);
	if (iReturn != -1)
	{
		return iReturn;
	}


	/*DUMP(gParams.mode);
	DUMP(gParams.files.size());
	DUMP(gParams.outPath);
	DUMP(gParams.calcMDefault);
	DUMP(gParams.calcMLumping);
	DUMP(gParams.solver);
	DUMP(gParams.quiet);
	DUMP(gParams.generatePNG);
	// ---
	DUMP(gParams.kappa);
	DUMP(gParams.k1);
	DUMP(gParams.k2);*/


	switch (gParams.mode)
	{
	case eHM_SIMPLE:
		iReturn = simple();	
		break;
	case eHM_TIMETEST:
		iReturn = timetest();
		break;
	case eHM_PRECISIONTEST:
		iReturn = precisiontest();
		break;
	case eHM_COMBINEDTEST:
		iReturn = timetest();
		if (iReturn != 0)
		{
			break;
		}
		iReturn = precisiontest();
		break;
	};

	return iReturn;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

int parseCmd(int argc, char* argv[])
{
	bool firstFile = true;

	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/helmholtz [-defM] [-lumpM] [-q] [-g] [-timetest] [-precisiontest] [-kappa " << gParams.kappa << "] [-k " << gParams.k1 << "," << gParams.k2  << "] [-o " << gParams.outPath << "] [" << gParams.files[0] << "] ..." << endl;
			cout << "  -defM           use mass lumping for the assembly of M" << endl;
			cout << "  -lumpM          use mass lumping for the assembly of M" << endl;
			cout << "  -q              don't run gnuplot to display plots" << endl;
			cout << "  -g              generate PNGs of the plots" << endl;
			cout << "  -timetest       measure and compare solvingtimes" << endl;
			cout << "  -precisiontest  measure and compare precisions" << endl;
			cout << "  -o              outpot folder for plots" << endl;
			cout << endl;
			cout << "If both '-defM' and '-lumpM' are used, comparisons are generated in the tests." << endl;
			cout << endl;
			cout << "Multiple files can be passed on UNIX systems using 'find', 'xargs' and pipelines:" << endl;
			cout << endl;
			cout << "  $ find data/mesh/carre?.msh | xargs bin/helmholtz -defM -lumpM -timetest" << endl;
			cout << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-defM") == 0)
		{
			gParams.calcMDefault = true;
		}
		else if (strcmp(argv[iArg], "-lumpM") == 0)
		{
			gParams.calcMLumping = true;
		}
		else if (strcmp(argv[iArg], "-q") == 0)
		{
			gParams.quiet = true;
		}
		else if (strcmp(argv[iArg], "-g") == 0)
		{
			gParams.generatePNG = true;
		}
		else if (strcmp(argv[iArg], "-timetest") == 0)
		{
			if (gParams.mode == eHM_PRECISIONTEST)
			{
				gParams.mode = eHM_COMBINEDTEST;
			}
			else
			{
				gParams.mode = eHM_TIMETEST;
			}
		}
		else if (strcmp(argv[iArg], "-precisiontest") == 0)
		{
			if (gParams.mode == eHM_TIMETEST)
			{
				gParams.mode = eHM_COMBINEDTEST;
			}
			else
			{
				gParams.mode = eHM_PRECISIONTEST;
			}
		}
		else if (strcmp(argv[iArg], "-o") == 0)
		{			
			iArg++;
			
			if (iArg < argc)
			{
				gParams.outPath = argv[iArg];
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
			if (firstFile)
			{
				gParams.files.clear();
				firstFile = false;
			}
			gParams.files.push_back(argv[iArg]);
		}
	}

	if (!gParams.calcMDefault && !gParams.calcMLumping)
	{
		// if no method was explicitly chosen, pick the default one
		gParams.calcMDefault = true;
	}

	return -1;
}

// ----------

bool loadMesh(const char* file, Mesh& mesh, clock_t& tLoadMesh)
{
	ifstream fileMesh(file);

	if (!fileMesh)
	{
		cout << "No valid mesh " << file << endl;
		return false;
	}

	RESETCLOCK();
    fileMesh >> mesh;
	CLOCK(tLoadMesh);

	return true;
}

void solveHelmholtz(Vector& u, const Mesh& mesh, bool lumping,
	clock_t& tMatA, clock_t& tMatM, clock_t& tMatB, clock_t& tRhsF, clock_t& tRhsG, clock_t& tSolve)
{
		// Assemble the matrices and vectors
		
		uint dim = mesh.countVertices();

		assert(dim == u.size());

		SparseMap Amap(dim, dim), Mmap(dim, dim), Bmap(dim, dim);
		Vector rhsF(dim), rhsG(dim);
		
		// A
		
		RESETCLOCK();
		Amap.constructA(mesh);
		CLOCK(tMatA);
		
		// M
		
		RESETCLOCK();
		if (lumping)
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
		
		switch(gParams.solver)
		{
		case eS_CG:
			RESETCLOCK();
			u = AMB.conjGradient(rhs);
			CLOCK(tSolve);
			break;
		case eS_JACOBI:
			RESETCLOCK();
			u = AMB.jacobi(rhs);
			CLOCK(tSolve);
			break;
		case eS_LU:
			// TODO: add LU solver to Sparse class if not there already
			// convert back to sparselil for use with lu
			/*RESETCLOCK();
			SparseLIL AMBlil(AMB,1);
			u = AMBlil.LU(rhs);
			CLOCK(tSolve);*/
			break;
		}

//	cout << "          Matrix NNZ: " << AMB.sizeNNZ() << " (" << (double)AMB.sizeNNZ() / (AMB.sizeRows() * AMB.sizeColumns()) * 100. << "%)" << endl;
}

int simple()
{
	cout << "Solving the problem over " << gParams.files.size() << " mesh(es)" << endl;

	// ---

	clock_t times[eT_END];

	RESETCLOCK();
	CLOCK(times[eT_START]);

	// Calculate approximate solutions and compare with the exact one

	for (uint iMsh = 0; iMsh < gParams.files.size(); iMsh++)
	{
		cout << endl << "Loading the mesh" << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		// ----------

		cout << "Calculating the exact solution..." << endl;

		Vector uExact(mesh.countVertices());
		uExact.constructFunc(mesh, helmholtz::u);

		{
			// visualization

			PlotMesh plotUExact("p1t2_uh_exacte", mesh, uExact, "Solution exacte");
			plotUExact.generate(ePT_MEDIT);
			plotUExact.generate(ePT_GNUPLOT_SURF, !gParams.quiet, false, "", "", 20);
			if (gParams.generatePNG) { plotUExact.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }
		}

		// ----------

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			cout << "Solving the linear system";
			if (lumping)
			{
				cout << " with lumping";
			}
			cout << "..." << endl;

			Vector u(mesh.countVertices());
			solveHelmholtz(u, mesh, lumping,
				times[eT_MATA], times[eT_MATM], times[eT_MATB], times[eT_RHSF], times[eT_RHSG], times[eT_SOLVE]);

			// ---
			
			// calculate errors

			Vector err = u - uExact;
		
			double error = err.norm2();
			double errorL2 = globalL2Error(mesh, uExact, u);
			double errorGradientL2 = globalL2GradError(mesh, uExact, u);

			// ----------

			printStats(gParams.files[iMsh].c_str(), mesh, times, error, errorL2, errorGradientL2);

			// ----------

			// visualization

			PlotMesh plotU("p1t2_uh_000", mesh, u, "Solution FEM");
			plotU.generate(ePT_MEDIT);
			plotU.generate(ePT_GNUPLOT_SURF, !gParams.quiet, false, "", "", 20);
			if (gParams.generatePNG) { plotU.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }
			
			PlotMesh plotErr("p1t2_err_000", mesh, err, "Erreur");
			plotErr.generate(ePT_MEDIT);
			plotErr.generate(ePT_GNUPLOT_SURF, !gParams.quiet, false, "", "", 20);
			if (gParams.generatePNG) { plotErr.generate(ePT_GNUPLOT_SURF, true, true, "", "", 20); }

			// ----------


			if (!lumping && gParams.calcMLumping)
			{
				// lumping is always the second step if we had the default method before
				lumping = true;
			}
			else
			{
				solved = true;
			}
		}
	}

	// ---

	CLOCK(times[eT_END]);

	return 0;
}

int timetest()
{
	cout << "Timetest over " << gParams.files.size() << " meshes" << endl;

	// ---

	clock_t times[eT_END];

	RESETCLOCK();
	CLOCK(times[eT_START]);

	// ---

	uint meshCount = gParams.files.size();

	Vector hMeshes(meshCount);

	Vector timesAssembleDefault(meshCount); 
	Vector timesSolveDefault(meshCount);

	Vector timesAssembleLumping(meshCount);
	Vector timesSolveLumping(meshCount);

	// Calculate solutions and clock the time

	for (uint iMsh = 0; iMsh < gParams.files.size(); iMsh++)
	{
		cout << endl << "Loading the mesh" << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		hMeshes(iMsh) = mesh.maxIncircleDiameter();

		// ----------

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			cout << "Solving the linear system";
			if (lumping)
			{
				cout << " with lumping";
			}
			cout << "..." << endl;

			Vector u(mesh.countVertices());
			solveHelmholtz(u, mesh, lumping,
				times[eT_MATA], times[eT_MATM], times[eT_MATB], times[eT_RHSF], times[eT_RHSG], times[eT_SOLVE]);

			// ----------

			//double assembleTime = ((double)(times[eT_MATA] + times[eT_MATM] + times[eT_MATB]
			//		+ times[eT_RHSF] + times[eT_RHSG]))	/ CLOCKS_PER_SEC;
			double assembleTime = ((double)times[eT_MATM]) / CLOCKS_PER_SEC;
			double solveTime = ((double) times[eT_SOLVE]) / CLOCKS_PER_SEC;

			if (lumping)
			{
				timesAssembleLumping(iMsh) = assembleTime;
				timesSolveLumping(iMsh) = solveTime;
			}
			else
			{
				timesAssembleDefault(iMsh) = assembleTime;
				timesSolveDefault(iMsh) = solveTime;
			}

			// ----------

			if (!lumping && gParams.calcMLumping)
			{
				// lumping is always the second step if we had the default method before
				lumping = true;
			}
			else
			{
				solved = true;
			}
		}
	}

	// ---

	// visualization

	if (gParams.calcMLumping && gParams.calcMDefault)
	{
		// compare default and lumping methods
		{
			Plot p("p2t4_timeSolve1", hMeshes, timesSolveDefault, "Time comparison (solving)", "",
					" w linespoints title 'default'");
			p.addYVector(timesSolveLumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Plot p("p2t4_timeAssemble1", hMeshes, timesAssembleDefault, "Time comparison (assembly)", "",
					" w linespoints title 'default'");
			p.addYVector(timesAssembleLumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}
	else
	{
		// plot default OR lumping times
		{
			Vector& y = gParams.calcMLumping ? timesSolveLumping : timesSolveDefault;
			string s = gParams.calcMLumping ? "Time with mass lumping (solving)" : "Time (solving)" ;

			Plot p("p2t4_timeSolve1", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector& y = gParams.calcMLumping ? timesAssembleLumping : timesAssembleDefault;
			string s = gParams.calcMLumping ? "Time with mass lumping (assembly)" : "Time (assembly)" ;

			Plot p("p2t4_timeSolve1", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}

	// ---

	CLOCK(times[eT_END]);

	return 0;
}

int precisiontest()
{
	cout << "Precisiontest over " << gParams.files.size() << " meshes" << endl;

	// ---

	clock_t times[eT_END];

	RESETCLOCK();
	CLOCK(times[eT_START]);

	// ---

	uint meshCount = gParams.files.size();

	Vector hMeshes(meshCount);

	Vector errorsDefault(meshCount);
	Vector errorsL2Default(meshCount);
	Vector errorsGradientL2Default(meshCount);

	Vector errorsLumping(meshCount);
	Vector errorsL2Lumping(meshCount);
	Vector errorsGradientL2Lumping(meshCount);

	// Calculate solutions and track the errors

	for (uint iMsh = 0; iMsh < gParams.files.size(); iMsh++)
	{
		cout << endl << "Loading the mesh" << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		hMeshes(iMsh) = mesh.maxIncircleDiameter();

		// ----------

		cout << "Calculating the exact solution..." << endl;

		Vector uExact(mesh.countVertices());
		uExact.constructFunc(mesh, helmholtz::u);

		// ----------

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			cout << "Solving the linear system";
			if (lumping)
			{
				cout << " with lumping";
			}
			cout << "..." << endl;

			Vector u(mesh.countVertices());
			solveHelmholtz(u, mesh, lumping,
				times[eT_MATA], times[eT_MATM], times[eT_MATB], times[eT_RHSF], times[eT_RHSG], times[eT_SOLVE]);

			// ----------

			// calculate errors

			Vector err = u - uExact;
		
			double error = err.norm2();
			double errorL2 = globalL2Error(mesh, uExact, u);
			double errorGradientL2 = globalL2GradError(mesh, uExact, u);

			if (lumping)
			{
				errorsLumping(iMsh) = error;
				errorsL2Lumping(iMsh) = errorL2;
				errorsGradientL2Lumping(iMsh) = errorGradientL2;
			}
			else
			{
				errorsDefault(iMsh) = error;
				errorsL2Default(iMsh) = errorL2;
				errorsGradientL2Default(iMsh) = errorGradientL2;
			}

			// ----------

			if (!lumping && gParams.calcMLumping)
			{
				// lumping is always the second step if we had the default method before
				lumping = true;
			}
			else
			{
				solved = true;
			}
		}
	}

	// ---

	// visualization

	if (gParams.calcMLumping && gParams.calcMDefault)
	{
		// compare default and lumping methods
		{
			Plot p("p2t4_error1", hMeshes, errorsDefault, "Error", "",
					" w linespoints title 'default'");
			p.addYVector(errorsLumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Plot p("p2t4_errorL21", hMeshes, errorsL2Default, "Error L2", "",
					" w linespoints title 'default'");
			p.addYVector(errorsL2Lumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Plot p("p2t4_errorGradientL21", hMeshes, errorsGradientL2Default, "Error Gradient L2", "",
					" w linespoints title 'default'");
			p.addYVector(errorsGradientL2Lumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}
	else
	{
		// plot default OR lumping precision
		{
			Vector& y = gParams.calcMLumping ? errorsLumping : errorsDefault;
			string s = gParams.calcMLumping ? "Errors with mass lumping" : "Errors" ;

			Plot p("p2t4_error1", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector& y = gParams.calcMLumping ? errorsL2Lumping : errorsL2Default;
			string s = gParams.calcMLumping ? "Errors L2 with mass lumping" : "Errors L2" ;

			Plot p("p2t4_errorL21", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector& y = gParams.calcMLumping ? errorsL2Lumping : errorsL2Default;
			string s = gParams.calcMLumping ? "Errors Gradient L2 with mass lumping" : "Errors Gradient L2" ;

			Plot p("p2t4_errorGradientL21", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "error");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}

	// ---

	CLOCK(times[eT_END]);
}

void printStats(const char* meshFile, const Mesh& mesh, clock_t times[],
				double error, double errorL2, double errorGradientL2)
{
	clock_t tCombinedAssembly = times[eT_MATA] + times[eT_MATM] + times[eT_MATB] + times[eT_RHSF] + times[eT_RHSG];

	cout << "---" << endl;
	
	LOGPARTTIME("A Assembly", times[eT_MATA], tCombinedAssembly);
	LOGPARTTIME("M Assembly", times[eT_MATM], tCombinedAssembly);
	LOGPARTTIME("B Assembly", times[eT_MATB], tCombinedAssembly);	
	LOGPARTTIME("F Assembly", times[eT_RHSF], tCombinedAssembly);
	LOGPARTTIME("G Assembly", times[eT_RHSG], tCombinedAssembly);

	cout << "---" << endl;
	
	LOGTIME("Combined Assembly", tCombinedAssembly);
	LOGTIME("Solving", times[eT_SOLVE]);
	LOGTIME("Load Mesh", times[eT_LOADMESH]);
		
	cout << "---" << endl;
	cout << "Computed: " << meshFile << " - Nv: " << mesh.countVertices() << 
	", Nt: " << mesh.countTriangles() << ", nE: " << mesh.countEdges() << ", h: " << mesh.maxIncircleDiameter() << endl;

	cout << "---" << endl;
	
	cout << "Error: " << error << endl;
	cout << "L2 error: " << errorL2 << endl;
	cout << "L2 grad error: " << errorGradientL2 << endl;
}

// ----------------------------------------------------------------------------

// FIN de la partie I :)
