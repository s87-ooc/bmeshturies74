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
//#include "mesh.h"
//#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	cout << "== Unittest ==" << endl;

	{
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
		cout << s << endl;

		Sparse t;

		t = s;
		cout << t << endl;
	}

	{
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
		cout << s << endl;

		// Jacobi test
		Vector b(2);
		b(0) = 11.;
		b(1) = 13.;

		// [ 2 1 ; 5 7 ] * [ x ; y ] = [ 11 ; 13 ]
		Vector solution = s.jacobi(b);
		cout << solution << endl;
	}

	{
		// [ 1 0 ; 0 1 ] 2x2

	    TVals vals;
		vals.push_back(1.);
		vals.push_back(1.);

	    TInd colInd;
		colInd.push_back(0);
		colInd.push_back(1);

		TInd rowPtr;
		rowPtr.push_back(0);
		rowPtr.push_back(1);
		rowPtr.push_back(2);

		Sparse s(vals, colInd, rowPtr, 2);
		cout << s << endl;

		// Jacobi test
		Vector b(2);
		b(0) = 11.;
		b(1) = 13.;

		// [ 1 0 ; 0 1 ] * [ x ; y ] = [ 11 ; 13 ]
		Vector solution = s.jacobi(b);
		cout << solution << endl;
	}

	{
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		cout << m << endl;
	}

	{
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

	return 0;
}
