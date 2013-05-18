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

/*

GENERAL STRUCTURE:

- function defs
- main with function calls
- function implementation

---

bin/wave data/mesh/cercle1.msh > plots (exact sol, position), times in console
	> simple
bin/wave -timetest
	> timetest
bin/wave -cfltest
	> precisiontest

---

Always do comparison over calculation methods

METHODS:	

bin/wave -lumpA
	> use lumping for A
bin/wave -defA
	> default construction for A

---

INPUT:
- 1 file
- 1 directory
- data/mesh/cercle*

>> Construction of filename list

---

OUTPUT:

bin/wave -o /dir
	> specify output directory
bin/wave -q
	> don't call gnuplot
bin/wave -g
	> generate graphics
bin/wave -r
	> generate plot at each step that can be combined into a video

---

PROBLEM SPECIFIC:

bin/wave -T
	> final time
bin/wave -dt
	> force time step
bin/wave -peak x,y
	> default position 0.1,0.1
	
*/

// ----------------------------------------------------------------------------

/** the program can loop over different meshes in different modes */
enum EWaveMode
{
	eWM_TIME = 0,	/** default mode, time is measured and dt picked optimally */
	eWM_CFL			/** we're varying dt to see stability */
};

// ----------------------------------------------------------------------------

struct SWaveParams
{
	EWaveMode mode;		/** mode for the main loop */
	double T;			/** maximal time */
	double dt;			/** time step */
	string fileMesh;	/** file of the current mesh */
	bool lumping;		/** enable mass lumping */
	bool render;		/** render a video of the time development */
	double framerate;	/** framerate of the video */
	uint meshCount;		/** number of meshes that should be tested automatically */
	uint dtCount;		/** how many dt values should be tested for CFL/Stability */
	bool save;			/** whether to save the plot */
	bool test;			/** select other set of files that were generated with freefem */
	// values below are used in test mode, [0] = normal [1] = mass lumping:
	double* hInner;			/** max inscribed diameter of the triangles */
	double* hOuter;			/** max circumscribed diameter of the triangles */
	double* dth;			/** dt/h for stability testing */
	double* times[2];		/** total computation time, [0] = normal [1] = mass lumping */
	Vector* uNorms[2];		/** norm of u after the final timestep, [0] = normal [1] = mass lumping */
	Vector* vNorms[2];		/** norm of v after the final timestep, [0] = normal [1] = mass lumping */
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
		return exp(-1. * (pow(v.x - 0.1, 2) + pow(v.y - 0.1, 2)) / 0.01);
	}
	
	/** initial value of the time derivative of u */
	double u1(const Vertex& v)
	{
		return 0.;
	}
};

// ----------------------------------------------------------------------------

void assembleIterationMatrices(const Mesh& mesh, Sparse& uMatU, Sparse& vMatU, Sparse& uMatV, Sparse& vMatV)
{
	uint dim = mesh.countVertices();

	// Assemble simple matrices

	SparseMap Amap(dim, dim), Mmap(dim, dim), Bmap(dim, dim);
	
	// A
	
	Amap.constructA(mesh);
	
	// M
	
	if (gParams.lumping)
	{
		Mmap.constructMlump(mesh);
	}
	else
	{
		Mmap.constructM(mesh);
	}
	
	// B
	
	Bmap.constructB(mesh);

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
	
	vMatV = Mmap;
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	clock_t t, tStart, tLoadMesh, tSolve, tEnd;
	clock_t* tSteps = 0;

	tStart = clock();

	// ---

	// commercial break
	
	cout << "  -- Wave (Partie 2), Projet Final [MM031], 2013" << endl;
	cout << "     Copyright (C) K. Podkanski, S. Stamenkovic" << endl;
	
	// ---
	
	// default values

	gParams.mode = eWM_TIME;
	gParams.T = 2.25;
	gParams.dt = -1.;
	gParams.fileMesh = "data/mesh/cercle1.msh";
	gParams.lumping = false;
	gParams.render = false;
	gParams.framerate = 25.;
	gParams.meshCount = 2;	// by default we'll treat data/mesh/cercle1.msh and data/mesh/cercle2.msh
	gParams.dtCount = 10;
	gParams.save = false;
	gParams.test = false;
	
	// get filenames and parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/wave [-lump] [-r fps] [-save] [-test] -T " << gParams.T << " -dt 0.001 "
			<< gParams.fileMesh << endl;
			cout << "       -lump = use mass lumping for the assembly of M" << endl;
			cout << "       -r fps = render a video of the time development (fps=25 if empty)" << endl;
			cout << "       -save = generate PNGs of the plots" << endl;
			cout << "       -test = perform tests on another set of meshes" << endl;
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
		else if (strcmp(argv[iArg], "-r") == 0)
		{
			gParams.render = true;
			
			if (iArg < argc - 1 && !strstr(argv[iArg+1], "-"))
			{
				iArg++;
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.framerate;
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
				cout << "Set dt manually to " << gParams.dt << endl;
			}
		}
		else
		{
			gParams.fileMesh = argv[iArg];
			gParams.dtCount = 1;
			gParams.meshCount = 1;
		}
	}

	// ----------
	
	gParams.hInner = new double[gParams.meshCount];
	gParams.hOuter = new double[gParams.meshCount];
	gParams.dth = new double[gParams.dtCount];
	gParams.times[0] = new double[gParams.meshCount];
	gParams.times[1] = new double[gParams.meshCount];
	gParams.uNorms[0] = new Vector[gParams.meshCount];
	gParams.uNorms[1] = new Vector[gParams.meshCount];
	gParams.vNorms[0] = new Vector[gParams.meshCount];
	gParams.vNorms[1] = new Vector[gParams.meshCount];
	
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	// fill the dt/h values, will be x-axis on CFL condition plots
	for (int i = 0; i < gParams.dtCount; i++)
	{
		gParams.dth[i] = 1.1 - 0.1 * i;
	}

	// iterate over the list of meshes
	
	for (int iMsh = 0; iMsh < gParams.meshCount; iMsh++)
	{
		// select the the default meshes or ffpp generated ones
	
		if (gParams.test)
		{
			// FFPP
			stringstream buf;
			buf << "data/mesh/cercle_ffpp_";
			buf << iMsh;
			buf >> gParams.fileMesh;
			
			gParams.fileMesh += ".msh";
		}
		else if (gParams.meshCount > 1)
		{
			stringstream buf;
			buf << "data/mesh/cercle";
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

		// calculate min/max diameters of the mesh triangles
		
		gParams.hInner[iMsh] = mesh.maxIncircleDiameter();
		gParams.hOuter[iMsh] = mesh.maxCircumcircleDiameter();

		// idx selects the [2] arrays to keep track of lumping values as well
		int idx = gParams.lumping ? 1 : 0;
		
		// we'll be starting in TIME mode, only one dt value is picked here
		int dtCount;
		
		if (gParams.mode == eWM_TIME)
		{
			dtCount = 1;
		}
		else
		{
			// prepare the vectors holding the norms for u and v
			gParams.uNorms[idx][iMsh] = Vector(gParams.dtCount);
			gParams.vNorms[idx][iMsh] = Vector(gParams.dtCount);
			dtCount = gParams.dtCount;
		}

		// iterate over dt when in CFL mode, else only one run
		
		for (uint iDt = 0; iDt < dtCount; iDt++)
		{
			// Prepare time steps
			
			if (gParams.mode == eWM_CFL)
			{
				// we're picking a fixed dth in CFT mode
				gParams.dt = gParams.dth[iDt] * gParams.hInner[iMsh];
			}
			else if (gParams.dt <= 0.)
			{
				gParams.dt = gParams.hInner[iMsh] / 3.;
			}
			
			uint nSteps = ceil(gParams.T / gParams.dt);
			tSteps = new clock_t[nSteps];

			// ----------
			
			// Assemble the matrices and vectors

			uint dim = mesh.countVertices();
			
			// vMatV = M
			Sparse uMatU, vMatU, uMatV, M;
			
			assembleIterationMatrices(mesh, uMatU, vMatU, uMatV, M);

			// ----------

			// solve the problem
			
			Vector u(dim);
			Vector uLast(dim);
			Vector v(dim);
			Vector vLast(dim);

			Vector originPos(nSteps+1);
			Vector timeVals(nSteps+1);
			
			// initial values
			uLast.constructFunc(mesh, wave::u0);
			vLast.constructFunc(mesh, wave::u1);

			// ---
			
			if (gParams.render)
			{
				PlotMesh plotD("wave_gnuplot_000", mesh, uLast);
				plotD.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl");
			}
			
			// ---
			
			originPos(0) = mesh.eval(0.0, 0.0, uLast);
			timeVals(0) = 0;
			
			cout << "Solving " << gParams.fileMesh << " with dt: " << gParams.dt << " T:" << " [0," << gParams.T << "]";
			if (gParams.lumping)
			{
				cout << " with mass lumping" << endl;
			}
			else
			{
				cout << endl;
			}
			
			// iteration over the time steps
			
			for (uint i = 0; i < nSteps; i++)
			{
				RESETCLOCK();

				M.newmark(u, v, uLast, vLast, uMatU, vMatU, uMatV);
				
				uLast = u;
				vLast = v;
				
				CLOCK(tSteps[i]);
				
				// in CFL test mode we may have divergency
				if (gParams.mode == eWM_CFL)
				{
					double _unorm = u.norm2();
					// the development after 2s should be way smaller than 10
					if (_unorm > 12.)
					{
						cout << "Divergency, abort" << endl;
						break;
					}
				}
				
				// track initial position
				originPos(i+1) = mesh.eval(0.0, 0.0, uLast);
				timeVals(i+1) = (i+1) * gParams.dt;
				
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
					
					PlotMesh plotD(plotFile.c_str(), mesh, u);
					plotD.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl");
				}
			}

			// computation time

			tSolve = 0;
			for (uint i = 0; i < nSteps; i++)
			{
				tSolve += tSteps[i];
			}
			
			double uNorm = u.norm2();
			double vNorm = v.norm2();
			
			if (gParams.mode == eWM_TIME)
			{
				gParams.times[idx][iMsh] = ((double)tSolve) / CLOCKS_PER_SEC;
			}
			else if (gParams.mode == eWM_CFL)
			{
				gParams.uNorms[idx][iMsh](iDt) = uNorm;
				gParams.vNorms[idx][iMsh](iDt) = vNorm;
			}
			
			// ---
			
			cout << "---" << endl;
			
			LOGTIME("Load Mesh", tLoadMesh);
			
			cout << "---" << endl;
			
			LOGTIME("Solving", tSolve);
			
			cout << "---" << endl;
			cout << "Computed: " << gParams.fileMesh << " - Nv: " << mesh.countVertices() << 
			", Nt: " << mesh.countTriangles() << ", nE: " << mesh.countEdges() << ", h: " << gParams.hInner[iMsh] << endl;

			cout << "---" << endl;
			
			cout << "norm2(u): " << uNorm << endl;
			cout << "norm2(v): " << vNorm << endl;
			cout << "dt/h: " << gParams.dt / gParams.hInner[iMsh] << endl;
			
			// ----------

			if (!gParams.test && !gParams.lumping && gParams.mode == eWM_TIME)
			{
				// Plot initial data if we're considering only a single mesh
				
				if (gParams.meshCount == 1)
				{
					PlotMesh plotU0("p2t1_u0", mesh, wave::u0, "Initial distribution");
					plotU0.generate(ePT_GNUPLOT, true);
					plotU0.generate(ePT_MEDIT);
					if (gParams.save) { plotU0.generate(ePT_GNUPLOT, true, true); }
					
					PlotMesh plotU1("p2t1_u1", mesh, wave::u1, "Initial velocity");
					plotU1.generate(ePT_GNUPLOT, true);
					plotU1.generate(ePT_MEDIT);
					if (gParams.save) { plotU1.generate(ePT_GNUPLOT, true, true); }
				}
				
				// plot last step
				
				ostringstream buf;
				buf << iMsh + 1;
				
				PlotMesh plotU(("p2t1_u_" + buf.str()).c_str(), mesh, u, "u(T, x)");
				plotU.generate(ePT_MEDIT);
				plotU.generate(ePT_GNUPLOT, true, false, "data/_gnuplot/wave.ptpl");
				if (gParams.save) { plotU.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl"); }
				
				PlotMesh plotV(("p2t1_v_" + buf.str()).c_str(), mesh, v, "v(T, x)");
				plotV.generate(ePT_MEDIT);
				plotV.generate(ePT_GNUPLOT, true, false, "data/_gnuplot/wave.ptpl");
				if (gParams.save) { plotV.generate(ePT_GNUPLOT, true, true, "data/_gnuplot/wave.ptpl"); }
				
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
				
				// Plot position of origin over time
				Plot plotOrigin(("p2t3_origin_" + buf.str()).c_str(), timeVals, originPos, "Position of origin over time", "", " w linespoints");
				plotOrigin.generate(ePT_GNUPLOT, true);
				if (gParams.save) { plotOrigin.generate(ePT_GNUPLOT, true, true); }
				
				// ---
				
				// render the time development as video if desired
				
				if (gParams.render)
				{
					PlotVideo::renderVideo("wave_gnuplot.avi", "wave_gnuplot_", gParams.framerate);
					break;
				}
			}
			
			// ----------
		
			SAFE_ARRDELETE(tSteps);
			
			// ----------
			
			if (gParams.mode == eWM_CFL)
			{
				cout << "CFL test: dt " << gParams.dt << "(" << iDt + 1 << "/" << gParams.dtCount << "), mesh "
					<< "(" << iMsh + 1 << "/" << gParams.meshCount << ")";
				if (gParams.lumping)
				{
					cout << " with mass lumping";
				}
				cout << endl;
			}
			
			if (gParams.render)
			{
				break;
			}
		}
		
		// ----------
			
		if (gParams.meshCount > 1 && iMsh == gParams.meshCount - 1)
		{
			if (!gParams.lumping)
			{
				// reset the iteration and calculate values with mass lumping
				cout << endl << "### Calculate values with mass lumping" << endl;
				iMsh = -1;
				gParams.lumping = true;
			}
			// switch the program mode
			else if (gParams.mode == eWM_TIME)
			{
				gParams.mode = eWM_CFL;
				gParams.lumping = false;
				cout << endl << endl << "=== CFL Test ===" << endl << endl;
				iMsh = -1;
			}
		}
		
		// reset the dt value, will be determined automatically on next step
		gParams.dt = -1.;
		
		if (gParams.render)
		{
			break;
		}
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
	
	// ----------
	
	// cleanup
	
	SAFE_ARRDELETE(gParams.hInner);
	SAFE_ARRDELETE(gParams.hOuter);
	SAFE_ARRDELETE(gParams.dth);
	SAFE_ARRDELETE(gParams.times[0]);
	SAFE_ARRDELETE(gParams.times[1]);
	SAFE_ARRDELETE(gParams.uNorms[0]);
	SAFE_ARRDELETE(gParams.uNorms[1]);
	SAFE_ARRDELETE(gParams.vNorms[0]);
	SAFE_ARRDELETE(gParams.vNorms[1]);

	tEnd = clock();
	
	cout << "---" << endl;
	LOGTIME("Total time", tEnd - tStart);
	
	return 0;
}

// FIN de la partie II :)
