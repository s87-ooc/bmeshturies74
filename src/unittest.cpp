/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
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
#include <math.h>

using namespace std;

// ----------------------------------------------------------------------------

// use this helper function only for very small sparse matrices!
template <class T>
void dump(const T& m)
{
	if (m.sizeColumns() > 10 || m.sizeRows() > 10)
	{
		cout << "Matrix is bigger than 10x10, use the << operator to print it" << endl;
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
		cout << "= Matrix LIL" << endl;
		
		SparseLIL m(3, 4);
	
		m(0, 1) = 2.0;
		m(0, 3) = 4.0;
		m(1, 1) = 1.0;
		m(2, 1) = 5.0;
		m(2, 0) = 3.0;
	
		dump(m);
		
		cout << "Rows: " << m.sizeRows() << endl;
		cout << "Columns: " << m.sizeColumns() << endl;
		cout << "NNZ: " << m.sizeNNZ() << endl;
		
		SparseLIL n(m);
		dump(n);
	}

	{
		cout << "= Matrix <-> Matrix LIL conversions" << endl;
		
		SparseLIL matLIL(3, 4);
	
		matLIL(0, 1) = 2.0;
		matLIL(0, 3) = 4.0;
		matLIL(1, 1) = 1.0;
		matLIL(2, 1) = 5.0;
		matLIL(2, 0) = 3.0;
		
		dump(matLIL);
		
		Sparse m(matLIL);
		dump(m);
		
		// ---
		
		TVals vals;
		vals.push_back(2.);
		vals.push_back(1.);
		vals.push_back(5.);
		vals.push_back(7.);

	    TInd colInd;
		colInd.push_back(0);
		colInd.push_back(1);
		colInd.push_back(0);
		colInd.push_back(2);

		TInd rowPtr;
		rowPtr.push_back(0);
		rowPtr.push_back(2);
		rowPtr.push_back(4);

		Sparse s(vals, colInd, rowPtr, 3);
		dump(s);
		
		SparseLIL nLIL(s);
		dump(nLIL);
	}
	
	
	{
		cout << "= Jacobi solver" << endl;
	
		// [ 2 1 ; 5 7 ] * [ x ; y ] = [ 11 ; 13 ]

		SparseLIL sLIL(2, 2);
		sLIL(0, 0) = 2.0;
		sLIL(0, 1) = 1.0;
		sLIL(1, 0) = 5.0;
		sLIL(1, 1) = 7.0;
		
		Sparse s(sLIL);
		dump(s);

		Vector b(2);
		b(0) = 11.;
		b(1) = 13.;

		Vector solution = s.jacobi(b);
		cout << solution << endl;
		
		cout << "Expected: [ 7.1111 ; -3.2222 ]" << endl;
	}

	{
		cout << "= Conjugate Gradient solver" << endl;
	
		{
			// [ 4 1 ; 1 3 ] * [ 0.0909 ; 0.6364 ] = [ 1 ; 2 ]
			
			SparseLIL sLIL(2, 2);
			sLIL(0, 0) = 4.0;
			sLIL(0, 1) = 1.0;
			sLIL(1, 0) = 1.0;
			sLIL(1, 1) = 3.0;
			
			Sparse s(sLIL);
			dump(s);

			Vector b(2);
			b(0) = 1.0;
			b(1) = 2.0;
			
			Vector solution = s.conjGradient(b);
			cout << solution << endl;
			
			cout << "Expected: [ 0.0909 ; 0.6364 ]" << endl;
		}
		
		{
			// [ 2 -1 0 ; -1 2 -1 ; 0 -1 2 ] * [ 2.5 ; 4 ; 3.5 ] = [ 1 ; 2 ; 3 ]
			
			SparseLIL sLIL(3, 3);
			sLIL(0, 0) = 2.0;
			sLIL(0, 1) = -1.0;
			sLIL(0, 2) = 0.0;
			sLIL(1, 0) = -1.0;
			sLIL(1, 1) = 2.0;
			sLIL(1, 2) = -1.0;
			sLIL(2, 0) = 0.0;
			sLIL(2, 1) = -1.0;
			sLIL(2, 2) = 2.0;
			
			Sparse s(sLIL);
			dump(s);

			Vector b(3);
			b(0) = 1.0;
			b(1) = 2.0;
			b(2) = 3.0;
			
			Vector solution = s.conjGradient(b);
			cout << solution << endl;
			
			cout << "Expected: [ 2.5 ; 4 ; 3.5 ]" << endl;
		}
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
		
        // ifstream meshfile;
        // meshfile.open("data/mesh/mesh.msh", ios::in);
        
        Mesh M("data/mesh/mesh.msh");
    }

    {
		cout << "= Mesh load" << endl;
		
        //Mesh mesh("data/mesh/mesh.msh");
		Mesh mesh("data/mesh/square_9.msh");

		/*for (uint v = 0; v < mesh.countVertices(); v++)
		{
			cout << v << ": " << mesh.V[v].x << " " << mesh.V[v].y << " (" << mesh.V[v].id << ")" << endl;
		}
		
		for (uint t = 0; t < mesh.countTriangles(); t++)
		{
			cout << t << ": " << mesh.T[t].V[0]->id << " " << mesh.T[t].V[1]->id << " " << mesh.T[t].V[2]->id << " (" << mesh.T[t].id << ")" << endl;
		}
		
		for (uint e = 0; e < mesh.countTriangles(); e++)
		{
			cout << e << ": " << mesh.E[e].V[0]->id << " " << mesh.E[e].V[1]->id << " (" << mesh.E[e].id << ")" << endl;
		}*/
		
		cout << "Mesh was loaded" << endl;
		
        uint Nv = mesh.countVertices();

        SparseMap A(Nv,Nv);
        SparseMap M(Nv,Nv);
        SparseMap B(Nv,Nv);

        A.constructA(mesh);
        
		// ---
		
		// A * (1, ..., 1) = 0
		
		Vector k(A.sizeRows());
		for (uint i = 0; i < k.size(); i++)
		{
			k(i) = 1.;
		}
		
		Sparse _A(A);
		Vector Ak = _A * k;

		cout << Ak << endl;
		
		for (uint i = 0; i < k.size(); i++)
		{
			for (uint j = 0; j < k.size(); j++)
			{
				cout << _A(i, j) << " ";
			}
			cout << endl;
		}

		cout << "Ak: " << Ak.norm2() << endl;
		
		return 0;
		
		// ---
		
		M.constructM(mesh);
        B.constructB(mesh);		
		
		// sum_i sum_j Mij = aire(omega)
		
		double area = 0.;
		for (uint t = 0; t < mesh.countTriangles(); t++)
		{
			area += mesh.T[t].area;
		}
		cout << "Area: " << area << endl;
		
		double areaM = 0.;
		for (uint i = 0; i < M.sizeRows(); i++)
		{
			for (uint j = 0; j < M.sizeColumns(); j++)
			{
				areaM += M(i,j);
			}
			//cout << i << endl;
		}
		cout << "Sum over Mij: " << areaM << endl;
		
		// ---
		
        double kappa = 10;
        SparseMap Res(Nv, Nv);

        M *= -pow(kappa,2);
        Res += A;
        Res += M;
        Res += B;

        Sparse R(Res);

        ofstream resultFile;
  		resultFile.open("data/result.dat");

  		resultFile << R;

  		resultFile.close();
    }
	
	{
		cout << "= Plot Mesh" << endl;
		
        Mesh M("data/mesh/square_9.msh");
		
		//PlotMesh p("plot01", "data/_gnuplot/surface.ptpl", &M);
		//p.generate(true);
    }

	return 0;
}
