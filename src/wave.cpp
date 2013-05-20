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

enum EWaveMode
{
	eWM_SIMPLE = 0,
	eWM_TIMETEST,
	eWM_CFLTEST,
	eWM_COMBINEDTEST
};

enum ETime
{
	eT_START = 0,
	eT_LOADMESH,
	eT_MATM,
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

/** simple run, calculate the solution over the specified mesh(es) and show position of origin */
int simple();

/** timetest, solve the problem on the meshes and compare calculation times, optionally with lumping */
int timetest();

/** cfltest, solve the problem over the meshes with varying timestep based on same dt/h factor */
int cfltest();

// ----------------------------------------------------------------------------

struct SWaveParams
{
	EWaveMode mode;
	std::vector<string> files;	/** list of mesh files we should perform our calculations on */
	string outPath;				/** path where generated plot files, etc should be stored in */
	bool calcMDefault;			/** whether to assemble A using the default method, can be combined with calcMLumping */
	bool calcMLumping;			/** whether to assemble A using mass lumping, can be combined with calcMDefault */
	ESolver solver;				/** which solver should be used for the linear system(s) */
	bool quiet;					/** if true, gnuplot won't be launched to display plots */
	bool generatePNG;			/** if true, gnuplot will be launched to generate a png */
	bool render;				/** render a video of the time development */
	double framerate;			/** framerate of the video */
	// ---
	double T;					/** maximal time */
	double dt;					/** time step, can be speficied  */
	std::vector<double> CFLs;	/** list of dt/h factors to run tests on */
	double peakX;				/** peak center of the initial function, x value */
	double peakY;				/** peak center of the initial function, y value */
};

/** global parameters */
SWaveParams gParams;

// ----------------------------------------------------------------------------

// we're keeping the problem specific functions in a namespace to avoid any ambiguity
namespace wave
{
	// NOTE: f and g are constantly 0 what simplifies the problem

	/** initial value of u */
	double u0(const Vertex& v)
	{
		return exp(-1. * (pow(v.x - gParams.peakX, 2) + pow(v.y - gParams.peakY, 2)) / 0.01);
	}

	/** initial value of the time derivative of u */
	double u1(const Vertex& v)
	{
		return 0.;
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

	cout << "  -- Wave (Partie 2), Projet Final [MM031], 2013" << endl;
	cout << "     Copyright (C) K. Podkanski, S. Stamenkovic" << endl << endl;

	// ---

	// default values

	gParams.mode = eWM_SIMPLE;
	gParams.files.push_back("data/mesh/cercle1.msh");
	gParams.outPath = "data/plots/";
	gParams.calcMDefault = false;	// will be set to true after parsing if no other method specified
	gParams.calcMLumping = false;
	gParams.solver = eS_CG;
	gParams.quiet = false;
	gParams.generatePNG = false;
	gParams.render = false;
	gParams.framerate = 25.;
	// ---
	gParams.T = 2.25;
	gParams.dt = -1.;	// forces automatic choice later
	gParams.CFLs.push_back(0.3);
	gParams.peakX = 0.1;
	gParams.peakY = 0.1;

	// ---

	// parse cmdline arguments, gives -1 if we should continue the program

	int iReturn = parseCmd(argc, argv);
	if (iReturn != -1)
	{
		return iReturn;
	}


	DUMP(gParams.mode);
	DUMP(gParams.files.size());
	DUMP(gParams.outPath);
	DUMP(gParams.calcMDefault);
	DUMP(gParams.calcMLumping);
	DUMP(gParams.solver);
	DUMP(gParams.quiet);
	DUMP(gParams.generatePNG);
	DUMP(gParams.render);
	DUMP(gParams.framerate);
	// ---
	DUMP(gParams.T);
	DUMP(gParams.dt);
	DUMP(gParams.CFLs.size());
	DUMP(gParams.peakX);
	DUMP(gParams.peakY);


	switch (gParams.mode)
	{
	case eWM_SIMPLE:
		iReturn = simple();
		break;
	case eWM_TIMETEST:
		iReturn = timetest();
		break;
	case eWM_CFLTEST:
		iReturn = cfltest();
		break;
	case eWM_COMBINEDTEST:
		iReturn = timetest();
		if (iReturn != 0)
		{
			break;
		}
		iReturn = cfltest();
		break;
	};

	return iReturn;
/*

				// plot last step

				ostringstream buf;
				buf << iMsh + 1;

				// Plot solving times if we're investigating a mesh in more detail
				if (gParams.meshCount == 1)
				{
					Vector x(nSteps);
					Vector y(nSteps);

					for (uint i = 0; i < nSteps; i++)
					{
						x(i) = i;
						y(i) = ((double)tSteps[i]) / CLOCKS_PER_SEC;
					}

					Plot p("p2t1_times", x, y, "solve time per step", "", " w linespoints");
					p.generate(ePT_GNUPLOT, true);
					if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
				}

				// ---

			}

	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	// plot combined data

	if (gParams.meshCount > 1)
	{
		{
			Vector x(gParams.dtCount);
			Vector y(gParams.dtCount);
			Vector ylump(gParams.dtCount);

			for (uint i = 0; i < gParams.dtCount; i++)
			{
				x(i) = gParams.dth[i];
				y(i) = gParams.uNorms[0][0](i);
			}

			Plot p("p2t2_cfl", x, y, "CFL condition", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector x(gParams.meshCount);
			Vector y(gParams.meshCount);
			Vector ylump(gParams.meshCount);

			for (uint i = 0; i < gParams.meshCount; i++)
			{
				x(i) = gParams.hInner[i];
				y(i) = gParams.times[0][i];
			}

			Plot p("p2t4_time2", x, y, "Mass Lumping: Time", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector x(gParams.dtCount);
			Vector y(gParams.dtCount);
			Vector ylump(gParams.dtCount);

			for (uint i = 0; i < gParams.dtCount; i++)
			{
				x(i) = gParams.dth[i];
				y(i) = gParams.uNorms[1][0](i);
			}

			Plot p("p2t4_stability", x, y, "Mass Lumping: Stability", "", " w linespoints");
			p.generate(ePT_GNUPLOT, true);
			if (gParams.save) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}
	*/
	return 0;
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
			cout << "Usage: bin/wave [-defM] [-lumpM] [-q] [-g] [-r] [-fps " << gParams.framerate << "] [-timetest] [-cfltest] [-cfls " << gParams.CFLs[0] << ",...] [-peak " << gParams.peakX << "," << gParams.peakY  << "] [-o " << gParams.outPath << "] [" << gParams.files[0] << "] ..." << endl;
			cout << "  -defM           use default method for the assembly of M" << endl;
			cout << "  -lumpM          use mass lumping for the assembly of M" << endl;
			cout << "  -q              don't run gnuplot to display plots" << endl;
			cout << "  -g              generate PNGs of the plots" << endl;
			cout << "  -r              render a video based on plots of each timestep" << endl;
			cout << "  -fps            framerate of the video to render" << endl;
			cout << "  -timetest       measure and compare solvingtimes" << endl;
			cout << "  -cfltest        measure the norm of the solution depending on cfl" << endl;
			cout << "  -cfls           list of cfl factors separated by comma" << endl;
			cout << "  -peak           position of the peak center for initial u" << endl;
			cout << "  -o              output folder for plots" << endl;
			cout << endl;
			cout << "If both '-defM' and '-lumpM' are used, comparisons are generated in the tests." << endl;
			cout << endl;
			cout << "Multiple files can be passed on UNIX systems using 'find', 'xargs' and pipelines:" << endl;
			cout << endl;
			cout << "  $ find data/mesh/cercle?.msh | xargs bin/wave -defM -lumpM -cfltest" << endl;
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
		else if (strcmp(argv[iArg], "-r") == 0)
		{
			gParams.render = true;
		}
		else if (strcmp(argv[iArg], "-fps") == 0)
		{
			iArg++;

			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.framerate;
			}
		}
		else if (strcmp(argv[iArg], "-timetest") == 0)
		{
			if (gParams.mode == eWM_CFLTEST)
			{
				gParams.mode = eWM_COMBINEDTEST;
			}
			else
			{
				gParams.mode = eWM_TIMETEST;
			}
		}
		else if (strcmp(argv[iArg], "-cfltest") == 0)
		{
			if (gParams.mode == eWM_TIMETEST)
			{
				gParams.mode = eWM_COMBINEDTEST;
			}
			else
			{
				gParams.mode = eWM_CFLTEST;
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
		else if (strcmp(argv[iArg], "-peak") == 0)
		{
			iArg++;

			if (iArg < argc)
			{
				char* sDim = strtok(argv[iArg], ",");

				{
					stringstream buf;
					buf << sDim;
					buf >> gParams.peakX;
				}

				sDim = strtok(0, ",");

				{
					stringstream buf;
					buf << sDim;
					buf >> gParams.peakY;
				}
			}
		}
		else if (strcmp(argv[iArg], "-cfls") == 0)
		{
			iArg++;

			if (iArg < argc)
			{
				char* sDim = strtok(argv[iArg], ",");

				while (sDim)
				{
					double cfl;
					stringstream buf;
					buf << sDim;
					buf >> cfl;

					gParams.CFLs.push_back(cfl);

					sDim = strtok(0, ",");
				}
			}
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

// ----------------------------------------------------------------------------

void printStats(const char* meshFile, const Mesh& mesh, clock_t times[],
				double uNorm, double vNorm, double timestep)
{
	double h = mesh.maxIncircleDiameter();

	cout << "---" << endl;

	LOGTIME("Assemble M", times[eT_MATM]);
	LOGTIME("Solving", times[eT_SOLVE]);
	LOGTIME("Load Mesh", times[eT_LOADMESH]);

	cout << "---" << endl;

	cout << "Computed: " << meshFile << " - Nv: " << mesh.countVertices() <<
	", Nt: " << mesh.countTriangles() << ", nE: " << mesh.countEdges() << ", h: " << h << endl;

	cout << "---" << endl;

	cout << "norm2(u): " << uNorm << endl;
	cout << "norm2(v): " << vNorm << endl;
	cout << "dt/h: " << timestep / h << endl;
}

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

void solveWave(Vector& uResult, Vector& vResult, Mesh& mesh, bool lumping,
	double timestep, bool origin, clock_t& tMatM, clock_t& tSolve)
{
	uint dim = mesh.countVertices();

	Vector u(dim);
	Vector v(dim);

	Vector uLast(dim);
	Vector vLast(dim);

	// ----------

	// vMatV = M
	Sparse uMatU, vMatU, uMatV, M;

	{
		// Assemble simple matrices

		SparseMap Amap(dim, dim), Mmap(dim, dim), Bmap(dim, dim);

		// A

		Amap.constructA(mesh);

		// M

		RESETCLOCK();
		if (gParams.calcMLumping)
		{
			Mmap.constructMlump(mesh);
		}
		else
		{
			Mmap.constructM(mesh);
		}
		CLOCK(tMatM);

		// B

		Bmap.constructB(mesh);

		// Assemble combined matrices

		SparseMap uMatUmap(dim, dim), vMatUmap(dim, dim), uMatVmap(dim, dim);

		SparseMap ABmap = Amap;
		ABmap += Bmap;

		uMatUmap = ABmap;
		uMatUmap *= -pow(timestep, 2) / 2.;
		uMatUmap += Mmap;
		uMatU = uMatUmap;

		vMatUmap = Mmap;
		vMatUmap *= timestep;
		vMatU = vMatUmap;

		uMatVmap = ABmap;
		uMatVmap *= -timestep / 2.;
		uMatV = uMatVmap;

		M = Mmap;
	}

	// ----------

	uint nSteps = ceil(gParams.T / timestep);

	clock_t* tSteps = new clock_t[nSteps];

	Vector originPos(nSteps+1);
	Vector timeVals(nSteps+1);

	// initial values
	uLast.constructFunc(mesh, wave::u0);
	vLast.constructFunc(mesh, wave::u1);

	if (gParams.render)
	{
		PlotMesh p("wave_gnuplot_000", mesh, uLast);
		p.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl");
	}

	// ---

	originPos(0) = mesh.eval(0.0, 0.0, uLast);
	timeVals(0) = 0;

	cout << "Solving with dt: " << timestep << " T:" << " [0," << gParams.T << "]";
	if (lumping)
	{
		cout << " with mass lumping";
	}
	cout << "..." << endl;

	// iteration over the time steps

	for (uint i = 0; i < nSteps; i++)
	{
		RESETCLOCK();

		M.newmark(u, v, uLast, vLast, uMatU, vMatU, uMatV);

		uLast = u;
		vLast = v;

		CLOCK(tSteps[i]);

		// in CFLtest mode we're handling divergency
		if (gParams.mode == eWM_CFLTEST)
		{
			double _unorm = u.norm2();
			// the development after 2.25s should not grow too much
			if (_unorm > 12.)
			{
				cout << "Divergency, abort" << endl;
				break;
			}
		}

		// track initial position
		originPos(i+1) = mesh.eval(0.0, 0.0, uLast);
		timeVals(i+1) = (i+1) * timestep;

		if (gParams.render)
		{
			stringstream buf;
			buf << "wave_gnuplot_";
			if (i+1 < 100)
			{
				if (i+1 < 10)
				{
					buf << "0";
				}
				buf << "0";
			}
			buf << i+1;

			string plotFile;
			buf >> plotFile;

			PlotMesh p(plotFile.c_str(), mesh, u);
			p.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl");
		}
	}

	// add up calculation time

	tSolve = 0;
	for (uint i = 0; i < nSteps; i++)
	{
		tSolve += tSteps[i];
	}

	// ---

	// Plot position of origin over time
	if (origin)
	{
		Plot p("p2t3_origin", timeVals, originPos, "Position of origin", "", " w linespoints");
		p.addScriptLine("set nokey");
		p.setAxisLabel(ePA_X, "t [s]");
		p.setAxisLabel(ePA_Y, "z");
		p.generate(ePT_GNUPLOT, !gParams.quiet);
		if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
	}

	// ---

	SAFE_ARRDELETE(tSteps);

	uResult = u;
	vResult = v;
}

// ----------------------------------------------------------------------------

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
		cout << endl << "Loading the mesh " << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		double h = mesh.maxIncircleDiameter();

		// ----------

		// setup timestep

		double dt = gParams.dt < 0 ? h * 0.3 : gParams.dt;

		// ----------

		// visualization of initial data

		{
			PlotMesh p("p2t1_u0", mesh, wave::u0, "Initial distribution");
			p.generate(ePT_GNUPLOT, true);
			p.generate(ePT_MEDIT);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			PlotMesh p("p2t1_u1", mesh, wave::u1, "Initial velocity");
			p.generate(ePT_GNUPLOT, true);
			p.generate(ePT_MEDIT);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		// ----------

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			Vector u(mesh.countVertices());
			Vector v(mesh.countVertices());

			solveWave(u, v, mesh, lumping, dt, true, times[eT_MATM], times[eT_SOLVE]);

			// ---

			// calculate norms

			double uNorm = u.norm2();
			double vNorm = v.norm2();

			// ----------

			printStats(gParams.files[iMsh].c_str(), mesh, times, uNorm, vNorm, dt);

			// ----------

			// visualization

			//plotU.generate(ePT_GNUPLOT_SURF, !gParams.quiet, false, "", "", 20);

			{
				PlotMesh p("p2t1_u", mesh, u, "u(T, x)");
				p.generate(ePT_MEDIT);
				p.generate(ePT_GNUPLOT, !gParams.quiet, false, "data/_gnuplot/wave.ptpl");
				if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl"); }
			}

			{
				PlotMesh p("p2t1_v", mesh, v, "v(T, x)");
				p.generate(ePT_MEDIT);
				p.generate(ePT_GNUPLOT, !gParams.quiet, false, "data/_gnuplot/wave.ptpl");
				if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl"); }
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


	// render the time development as video if desired

	if (gParams.render)
	{
		PlotVideo::renderVideo("wave_gnuplot.avi", "wave_gnuplot_", gParams.framerate);
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
		cout << endl << "Loading the mesh " << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		hMeshes(iMsh) = mesh.maxIncircleDiameter();

		// ----------

		// setup timestep

		double dt = gParams.dt < 0 ? hMeshes(iMsh) * 0.3 : gParams.dt;

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			Vector u(mesh.countVertices());
			Vector v(mesh.countVertices());

			solveWave(u, v, mesh, lumping, dt, false, times[eT_MATM], times[eT_SOLVE]);

			// ----------

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
			Plot p("p2t4_timeSolve2", hMeshes, timesSolveDefault, "Time comparison (solving)", "",
					" w linespoints title 'default'");
			p.addYVector(timesSolveLumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Plot p("p2t4_timeAssemble2", hMeshes, timesAssembleDefault, "Time comparison (assembly)", "",
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

			Plot p("p2t4_timeSolve2", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector& y = gParams.calcMLumping ? timesAssembleLumping : timesAssembleDefault;
			string s = gParams.calcMLumping ? "Time with mass lumping (assembly)" : "Time (assembly)" ;

			Plot p("p2t4_timeSolve2", hMeshes, y, s.c_str(), "", " w linespoints");
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

int cfltest()
{
	cout << "CFLtest over " << gParams.files.size() << " meshes and " << gParams.CFLs.size() << " factors" << endl;
/*
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
		cout << endl << "Loading the mesh " << gParams.files[iMsh] << "..." << endl;

		Mesh mesh;
		if (!loadMesh(gParams.files[iMsh].c_str(), mesh, times[eT_LOADMESH]))
		{
			return 1;
		}

		hMeshes(iMsh) = mesh.maxIncircleDiameter();

		// ----------

		// setup timestep

		double dt = gParams.dt < 0 ? hMeshes(iMsh) / .3 : gParams.dt;

		// problem will be solved with and/or without mass lumping (method), no lumping first if desired

		bool solved = false;
		bool lumping = !gParams.calcMDefault;

		while (!solved)
		{
			Vector u(mesh.countVertices());
			Vector v(mesh.countVertices());

			solveWave(u, v, mesh, lumping, dt, false, times[eT_MATM], times[eT_SOLVE]);

			// ----------

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
			Plot p("p2t4_timeSolve2", hMeshes, timesSolveDefault, "Time comparison (solving)", "",
					" w linespoints title 'default'");
			p.addYVector(timesSolveLumping, " w linespoints title 'mass lumping'");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Plot p("p2t4_timeAssemble2", hMeshes, timesAssembleDefault, "Time comparison (assembly)", "",
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

			Plot p("p2t4_timeSolve2", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}

		{
			Vector& y = gParams.calcMLumping ? timesAssembleLumping : timesAssembleDefault;
			string s = gParams.calcMLumping ? "Time with mass lumping (assembly)" : "Time (assembly)" ;

			Plot p("p2t4_timeSolve2", hMeshes, y, s.c_str(), "", " w linespoints");
			p.setAxisLabel(ePA_X, "h");
			p.setAxisLabel(ePA_Y, "t [s]");
			p.addScriptLine("set nokey");
			p.generate(ePT_GNUPLOT, !gParams.quiet);
			if (gParams.generatePNG) { p.generate(ePT_GNUPLOT, true, true); }
		}
	}

	// ---

	CLOCK(times[eT_END]);
*/
	return 0;
}

// FIN de la partie II :)
