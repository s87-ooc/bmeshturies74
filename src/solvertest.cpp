/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: solvertest.cpp

 Description: mesure du temps des solveurs et génération des matrices

**************************************************************/

#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iomanip>

// ---

#include "algebra.h"
#include "visualization.h"
#include "clock.h"

using namespace std;

// ----------------------------------------------------------------------------

struct SGenerationParams
{
	bool regenerate;	/** whether example matrices should be regenerated */
	double density;		/** nnz / (rows * columns), ratio the same for all matrices, minimum not singular */
	uint count;			/** how many example matrices we have */
	uint* dimensions;	/** dynamic array of example matrix dimensions */
	double* times;		/** dynamic array of computation times, 3*dimension: conjGradient, jacobi, LU */
	// ---
	Sparse* matrices;	/** dynamic array of matrices of the linear system */
	Vector* rhs;		/** dynamic array of rhs vectors of the linear system */
};

/** global parameters */
SGenerationParams gParams;

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	// ---
	
	// default values

	gParams.regenerate = false;
	gParams.density = 0.001;	// similar to problem matrices of carre2.msh/carre3.msh/carre4.msh
	gParams.count = 4;
	gParams.dimensions = 0;
	gParams.times = 0;
	gParams.matrices = 0;
	gParams.rhs = 0;
	
	// get parameters from cmdline arguments (if any)
	
	for (int iArg = 1; iArg < argc; iArg++)
	{
		if (strcmp(argv[iArg], "-h") == 0)
		{
			cout << "Usage: bin/solvertest [-r] -density " << gParams.density << " -count " << gParams.count << " dim1,...,dim[count]" << endl;
			cout << "       -r = regenerate existing files" << endl;
			return 0;
		}
		else if (strcmp(argv[iArg], "-r") == 0)
		{
			gParams.regenerate = true;
		}
		else if (strcmp(argv[iArg], "-density") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.density;
			}
		}
		else if (strcmp(argv[iArg], "-count") == 0)
		{
			iArg++;
			
			if (iArg < argc)
			{
				stringstream buf;
				buf << argv[iArg];
				buf >> gParams.count;
			}
		}
		else
		{
			if (iArg < argc)
			{
				gParams.dimensions = new uint[gParams.count];
				
				char* dim = strtok(argv[iArg], ",");
				for (uint i = 0; i < gParams.count; i++)
				{
					assert(dim);
				
					stringstream buf;
					buf << dim;
					buf >> gParams.dimensions[i];
				
					dim = strtok(0, ",");
				}
			}
		}
	}

	cout << "== Partie I: Solvertest ==" << endl;
	
	// ----------
	
	if (!gParams.dimensions)
	{
		gParams.dimensions = new uint[gParams.count];
		for (uint i = 0; i < gParams.count; i++)
		{
			gParams.dimensions[i] = ceil(pow(10, i+1));
		}
	}
	
	cout << "Dimensions: ";
	DUMP_ARR(gParams.dimensions, gParams.count);
	
	gParams.times = new double[gParams.count * 3];
	gParams.matrices = new Sparse[gParams.count];
	gParams.rhs = new Vector[gParams.count];
	
	// ----------
	
	// generate linear systems if needed
	
	// prepare the random number generator
	
	for (uint k = 0; k < gParams.count; k++)
	{
		// generate positive definite matrices over the product of a lower triangular matrix and its transpose
		// fill the diagonal with positive values to assure that its invertible
		
		SparseLIL m(gParams.dimensions[k], gParams.dimensions[k]);
		for (uint i = 0; i < m.sizeRows(); i++)
		{
			for (uint j = 0; j < m.sizeColumns(); j++)
			{
				
				
				// make sure we have a positive value on the diagonal
				while (m(j, j) < EQ_TOL)
				{
					m(j, j) = j + 1;
				}
			}
		}
		
		gParams.matrices[k] = m.prodTranspose();

		// generate solution vectors
	
		Vector x = Vector(gParams.dimensions[k]);
		gParams.rhs[k] = gParams.matrices[k] * x;
	}
	
	// ----------
	
	// solve systems and clock the times
	
	clock_t tSolve;
	for (uint k = 0; k < gParams.count; k++)
	{
		// Conjugate Gradient
	
		RESETCLOCK();
		
		//gParams.matrices[k].conjGradient(gParams.rhs[k]);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3] = ((double)tSolve) / CLOCKS_PER_SEC;
		
		// Jacobi
		
		RESETCLOCK();
		
		//gParams.matrices[k].jacobi(gParams.rhs[k]);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3 + 1] = ((double)tSolve) / CLOCKS_PER_SEC;
		
		// LU
		
		RESETCLOCK();
		
		//gParams.matrices[k].lu(gParams.rhs[k]);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3 + 2] = ((double)tSolve) / CLOCKS_PER_SEC;
	}
	
	// ----------

	// write logfile
	
	{
		ofstream fileClock("data/linsys/times_conjGradient.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3] << endl;
		}
	}

	{
		ofstream fileClock("data/linsys/times_jacobi.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3 + 1] << endl;
		}
	}
	
	{
		ofstream fileClock("data/linsys/times_lu.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3 + 2] << endl;
		}
	}
	
	// ----------
	
	// cleanup
	
	SAFE_ARRDELETE(gParams.dimensions);
	SAFE_ARRDELETE(gParams.times);
	SAFE_ARRDELETE(gParams.matrices);
	SAFE_ARRDELETE(gParams.rhs);

	return 0;
}