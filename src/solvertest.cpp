/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenković (stjepan@stjepan.net)

 ---

     Fichier: solvertest.cpp

 Description: mesure du temps des solveurs et génération des matrices

**************************************************************/

#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <cstdlib>

// ---

#include "algebra.h"
#include "visualization.h"
#include "clock.h"

using namespace std;

// ----------------------------------------------------------------------------

struct SGenerationParams
{
	bool regenerate;	/** whether example linear systems should be regenerated */
	uint count;			/** how many example systems we have */
	double density;		/** nnz / (rows * columns), ratio the same for all matrices, plus diagonal nnz */
	double minEntry;	/** minimal entry of randomly generated values */
	double maxEntry;	/** maximal entry of randomly generated values */
	uint* dimensions;	/** dynamic array of test matrix dimensions */
	double* times;		/** dynamic array of computation times, 3 * count: conjGradient, jacobi, LU */
	double* errors;		/** dynamic array of computation errors, 3 * count: conjGradient, jacobi, LU */
};

/** global parameters */
SGenerationParams gParams;

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	// default values

	gParams.regenerate = false;
	gParams.count = 2;
	gParams.density = 0.001;	// similar to problem matrices of carre2.msh/carre3.msh/carre4.msh
	gParams.minEntry = -3.;
	gParams.maxEntry = 8.;
	gParams.dimensions = 0;
	gParams.times = 0;
	gParams.errors = 0;
	
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
			// dim1,...,dim[count]
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

	cout << "== Partie I: Solvertest ==" << endl;
	
	// ----------
	
	if (!gParams.dimensions)
	{
		gParams.dimensions = new uint[gParams.count];
		for (uint i = 0; i < gParams.count; i++)
		{
			// default dimensions are 10^i
			gParams.dimensions[i] = ceil(pow(10, i+1));
		}
	}
	
	cout << "Dimensions: ";
	DUMP_ARR(gParams.dimensions, gParams.count);
	
	gParams.times = new double[gParams.count * 3];
	gParams.errors = new double[gParams.count * 3];
	
	// ----------

	// prepare the random number generator	
	// TODO: use C++11 random generators with actual uniform distribution
	srand(time(0));
	double rangeEntry = gParams.maxEntry - gParams.minEntry;
	int NZ_MAX = (int)(gParams.density * (double)RAND_MAX);
	double ratio = rangeEntry / RAND_MAX;
	
	// ----------
	
	// iterate over systems we need to calculate
	
	for (uint k = 0; k < gParams.count; k++)
	{
		uint dim = gParams.dimensions[k];
		Sparse A;
		Vector rhs(dim);
		Vector solution(dim);
		Vector exact_solution(dim);
		
		// construct the filename of the linear system

		stringstream buf;
		buf << "data/linsys/test_" << gParams.dimensions[k] << ".linsys";
		string fn;
		buf >> fn;

		ifstream fileLinSys(fn.c_str());
		
		// (re)generate linear system if needed
		
		if (gParams.regenerate || !fileLinSys)
		{
			if (fileLinSys)
			{
				fileLinSys.close();
			}
		
			// generate positive definite matrices over the product of a lower triangular matrix
			// and its transpose, fill the diagonal with positive values to assure that its invertible
			
			SparseLIL m(dim, dim);

			for (uint i = 0; i < m.sizeRows(); i++)
			{
				for (uint j = 0; j < i; j++)
				{
					if (rand() < NZ_MAX)
					{
						m(i, j) = rand() * ratio + gParams.minEntry;
					}
				}
				// make sure we have a positive value on the diagonal
				while (m(i, i) < EQ_TOL)
				{
					m(i, i) = rand() * ratio;
				}
				
				if (i % 100 == 0)
					cout << i << endl;
			}
			cout << dim << ": matrice triangulaire inférieure generée" << endl;
			
			A = m.prodTranspose();

			cout << dim << ": matrice generée" << endl;
			
			// generate solution vectors
		
			exact_solution = Vector(dim);
			for (uint i = 0; i < dim; i++)
			{
				exact_solution(i) = rand() * ratio + gParams.minEntry;
			}
			
			rhs = A * exact_solution;
			
			ofstream f(fn.c_str());
			f << A << rhs << exact_solution;
		}
		else
		{
			fileLinSys >> A >> rhs >> exact_solution;
			fileLinSys.close();
		}
	
		cout << dim << " preparée" << endl;
	
		// ----------
		
		// solve systems and clock the times
		
		clock_t tSolve;
		
		// Conjugate Gradient
	
		RESETCLOCK();
		
		solution = A.conjGradient(rhs);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3] = ((double)tSolve) / CLOCKS_PER_SEC;
		gParams.errors[k * 3] = (solution - exact_solution).norm2();
		
		cout << dim << " Gradient Conjugué: " << gParams.times[k * 3] << "s " << gParams.errors[k * 3] << endl;
		
		// Jacobi
		
		RESETCLOCK();
		
		solution = A.jacobi(rhs);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3 + 1] = ((double)tSolve) / CLOCKS_PER_SEC;
		gParams.errors[k * 3 + 1] = (solution - exact_solution).norm2();
		
		cout << dim << " Jacobi: " << gParams.times[k * 3 + 1] << "s " << gParams.errors[k * 3 + 1] << endl;
		
		// LU
		
		RESETCLOCK();
		
		solution = A.LU(rhs);
		
		CLOCK(tSolve);
		
		gParams.times[k * 3 + 2] = ((double)tSolve) / CLOCKS_PER_SEC;
		gParams.errors[k * 3 + 2] = (solution - exact_solution).norm2();
		
		cout << dim << " LU: " << gParams.times[k * 3 + 2] << "s " << gParams.errors[k * 3 + 2] << endl;
	}
	
	// ----------

	// write logfiles
	
	{
		ofstream fileClock("data/linsys/times_conjGradient.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3] << endl;
		}
		
		ofstream fileError("data/linsys/errors_conjGradient.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileError << gParams.dimensions[k] << " " << gParams.errors[k * 3] << endl;
		}
	}

	{
		ofstream fileClock("data/linsys/times_jacobi.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3 + 1] << endl;
		}
		
		ofstream fileError("data/linsys/errors_jacobi.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileError << gParams.dimensions[k] << " " << gParams.errors[k * 3 + 1] << endl;
		}
	}
	
	{
		ofstream fileClock("data/linsys/times_lu.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileClock << gParams.dimensions[k] << " " << gParams.times[k * 3 + 2] << endl;
		}
		
		ofstream fileError("data/linsys/errors_lu.log");
		for (uint k = 0; k < gParams.count; k++)
		{
			fileError << gParams.dimensions[k] << " " << gParams.errors[k * 3 + 2] << endl;
		}
	}
	
	// ----------

	// plot results
	
	// ----------
	
	// cleanup
	
	SAFE_ARRDELETE(gParams.dimensions);
	SAFE_ARRDELETE(gParams.times);
	SAFE_ARRDELETE(gParams.errors);

	return 0;
}