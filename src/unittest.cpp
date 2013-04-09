/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: unittest.cpp

 Description: test

**************************************************************/

#include <iostream>
#include <fstream>

#include "algebra.h"
#include "mesh.h"
#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

// use this helper function only for very small sparse matrices!
void dump(Sparse& m)
{
	if (m.sizeColumns() > 10 || m.sizeRows() > 10)
	{
		cout << "Matrix is bigger than 9x9, use the << operator to print it" << endl;
	}
	
	cout << "[ ";
	for (uint i = 0; i < m.sizeRows(); i++)
	{
		for (uint j = 0; j < m.sizeColumns(); j++)
		{
			cout << m(i, j) << " ";
		}
		
		if (i < m.sizeRows() - 1)
		{
			cout << "; ";
		}
	}
	cout << "] " << m.sizeRows() << "x" << m.sizeColumns() << endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	cout << "== Unittest ==" << endl;

	{
		cout << "= Vector" << endl;
		
		Vector v(3);
		
		cout << v << endl;

		v(0) = 1.0;
		v(1) = 2.0;
		v(2) = 3.0;

		cout << v << endl;

		Vector w;

		w = v;
		cout << w << endl;

		cout << v.norm2() << endl;

		w = 3.0 * v;
		cout << w << endl;

		w = v * 3.0;
		cout << w << endl;
	}

	{
		cout << "= Matrix construcor" << endl;
		
		// [ 2 1 ; 5 7 ] 2x2

	    TVals vals;
		vals.push_back(2.);
		vals.push_back(1.);
		vals.push_back(5.);
		vals.push_back(7.);

	    TInd colInd;
		colInd.push_back(0);
		colInd.push_back(1);
		colInd.push_back(0);
		colInd.push_back(1);

		TInd rowPtr;
		rowPtr.push_back(0);
		rowPtr.push_back(2);
		rowPtr.push_back(4);

		Sparse s(vals, colInd, rowPtr, 2);
		dump(s);

		Sparse t;

		t = s;
		dump(t);
	}

	{
		cout << "= Jacobi solver" << endl;
	
		// [ 2 1 ; 5 7 ] * [ x ; y ] = [ 11 ; 13 ]

	    TVals vals;
		vals.push_back(2.);
		vals.push_back(1.);
		vals.push_back(5.);
		vals.push_back(7.);

	    TInd colInd;
		colInd.push_back(0);
		colInd.push_back(1);
		colInd.push_back(0);
		colInd.push_back(1);

		TInd rowPtr;
		rowPtr.push_back(0);
		rowPtr.push_back(2);
		rowPtr.push_back(4);

		Sparse s(vals, colInd, rowPtr, 2);
		//dump(s);

		Vector b(2);
		b(0) = 11.;
		b(1) = 13.;

		Vector solution = s.jacobi(b);
		cout << solution << endl;
		
		cout << "Expected: [ 7.1111 ; -3.2222 ]" << endl;
	}

	{
		cout << "= Conjugate Gradient solver" << endl;
	
		// [ 4 1 ; 1 3 ] * [ 0.0909 ; 0.6364 ] = [ 1 ; 2 ]
		
	    TVals vals;
		vals.push_back(4.);
		vals.push_back(1.);
		vals.push_back(1.);
		vals.push_back(3.);

	    TInd colInd;
		colInd.push_back(0);
		colInd.push_back(1);
		colInd.push_back(0);
		colInd.push_back(1);

		TInd rowPtr;
		rowPtr.push_back(0);
		rowPtr.push_back(2);
		rowPtr.push_back(4);

		Sparse s(vals, colInd, rowPtr, 2);
		dump(s);

		Vector b(2);
		b(0) = 1.0;
		b(1) = 2.0;
		
		Vector solution = s.conjGradient(b);
		cout << solution << endl;
		
		cout << "Expected: [ 0.0909 ; 0.6364 ]" << endl;
	}

	{
		cout << "= Matrix load" << endl;
	
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		dump(m);
	}

	{
		cout << "= Matrix multiplication" << endl;
	
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		Vector v(3);
		v(0) = 1.;
		v(1) = 2.;
		v(2) = 3.;
		
		Vector w = m * v;
		
		cout << w << endl;
	}

	{
		cout << "= Matrix product v' * M * w" << endl;
		
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		Vector v(3);
		v(0) = 1.;
		v(1) = 2.;
		v(2) = 3.;

		Vector w(3);
		w(0) = 3.;
		w(1) = 4.;
		w(2) = 1.;

		double d = prod(v, m, w);
		
		cout << d << endl;
	}
    
    {
		cout << "= Mesh load" << endl;
		
        ifstream meshfile;
        meshfile.open("data/mesh/mesh.msh", ios::in);
        
        Mesh M;
        meshfile >> M;
    }
	
	{
		cout << "= Plot Mesh" << endl;
		
        ifstream meshfile;
        meshfile.open("data/mesh/mesh.msh", ios::in);
        
        Mesh M;
        meshfile >> M;
		
		Plot p("plot01", "data/_gnuplot/surface.ptpl", &M);
		//p.generate(true);
    }

	return 0;
}
