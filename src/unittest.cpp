/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenković (stjepan@stjepan.net)

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

#define TESTTITLE(T) cout << endl << "= " << T << " =" << endl << endl

template<class T>
void _writeMatrix(const char* fn, const T& m)
{
	ofstream os(fn);
	for (uint i = 0; i < m.sizeRows(); i++)
	{
		for (uint j = 0; j < m.sizeColumns(); j++)
		{
			os << m(i, j) << endl;
		}
	}
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	cout << "== Unittest ==" << endl;


	{
		TESTTITLE("Matrix assembly");
        
        Mesh mesh("data/mesh/square_9.msh");
		
		for (uint i = 0; i < mesh.countVertices(); i++)
		{
			cout << i+1 << " V) " 
				<< mesh.V[i].x << " " << mesh.V[i].y << " " /*<< mesh.V[i].id + 1*/ << endl;
		}
		
		for (uint i = 0; i < mesh.countTriangles(); i++)
		{
			cout << i+1 << " T) " 
				<< mesh.T[i](0).id+1 << "-" << mesh.T[i](1).id+1 << "-" << mesh.T[i](2).id+1
				<< " " << mesh.T[i].area << endl;
		}
		
		for (uint i = 0; i < mesh.countEdges(); i++)
		{
			cout << i+1 << " E) " 
				<< mesh.E[i](0).id+1 << "-" << mesh.E[i](1).id+1
				<< " " << mesh.E[i].length << endl;
		}
		
		uint dim = mesh.countVertices();
		
		SparseMap Amap(dim, dim), Bmap(dim, dim), Mmap(dim, dim);
		
		Amap.constructA(mesh);
		DUMP_MAT(Amap);
		
		Bmap.constructB(mesh);
		DUMP_MAT(Bmap);
		
		Mmap.constructM(mesh);
		DUMP_MAT(Mmap);
		
		// we can load this in Scilab with "A = read('A.dat',9,9)"
		//_writeMatrix("A.dat", Amap);
		//_writeMatrix("B.dat", Bmap);
		//_writeMatrix("M.dat", Mmap);
    }
	
	{
		TESTTITLE("Vector");
		
		Vector v(3);
		
		DUMP_VEC(v);

		v(0) = 1.0;
		v(1) = 2.0;
		v(2) = 3.0;

		DUMP_VEC(v);

		Vector w;

		w = v;
		DUMP_VEC(w);

		cout << v.norm2() << endl;

		w = 3.0 * v;
		DUMP_VEC(w);

		w = v * 3.0;
		DUMP_VEC(w);
	}
	
	{
		TESTTITLE("Vector I/O");
		
		Vector v(3);
		v(0) = 2.;
		v(1) = -1.;
		v(2) = 42.;
		
		Vector w(2);
		w(0) = 9.;
		w(1) = 7.;
		
		string fn = "_vectest.dat";
		
		{
			ofstream fileOut(fn.c_str());
			fileOut << v << w;
		}
		
		Vector x, y;
		
		{
			ifstream fileIn(fn.c_str());
			fileIn >> x >> y;			
		}
		
		DUMP_VEC(x);
		DUMP_VEC(y);
		
		remove(fn.c_str());
	}

	{
		TESTTITLE("Sparse");
		
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
		DUMP_MAT(s);

		Sparse t;

		t = s;
		DUMP_MAT(t);
	}
	
	{
		TESTTITLE("Sparse I/O");
		
		SparseLIL _m(3, 3);
		_m(0, 0) = 1.;
		_m(1, 0) = 4.;
		_m(1, 1) = 2.;
		_m(2, 0) = 5.;
		_m(2, 2) = 3.;

		SparseLIL _n(3, 4);
		_n(0, 1) = 2.0;
		_n(0, 3) = 4.0;
		_n(1, 1) = 1.0;
		_n(2, 1) = 5.0;
		_n(2, 0) = 3.0;
		
		Sparse m(_m);
		Sparse n(_n);

		DUMP_MAT(m);
		DUMP_MAT(n);
		
		string fn = "_mattest.dat";
		
		{
			ofstream fileOut(fn.c_str());
			fileOut << m << n;
		}
		
		Sparse a, b;
		
		{
			ifstream fileIn(fn.c_str());
			fileIn >> a >> b;			
		}
		
		DUMP_MAT(a);
		DUMP_MAT(b);
		
		remove(fn.c_str());
	}

	{
		TESTTITLE("Sparse LIL");
		
		SparseLIL m(3, 4);
	
		m(0, 1) = 2.0;
		m(0, 3) = 4.0;
		m(1, 1) = 1.0;
		m(2, 1) = 5.0;
		m(2, 0) = 3.0;
	
		DUMP_MAT(m);
		
		cout << "Rows: " << m.sizeRows() << endl;
		cout << "Columns: " << m.sizeColumns() << endl;
		cout << "NNZ: " << m.sizeNNZ() << endl;
		
		SparseLIL n(m);
		DUMP_MAT(n);
	}

	{
		TESTTITLE("matrix * transpose(matrix)");
		
		SparseLIL m(3, 3);
	
		m(0, 0) = 1.;
		m(1, 0) = 4.;
		m(1, 1) = 2.;
		m(2, 0) = 5.;
		m(2, 2) = 3.;

		DUMP_MAT(m);

		SparseLIL mprod = m.prodTranspose();
		
		DUMP_MAT(mprod);
		
		cout << "Expected: [ 1 4 5 ; 4 20 20 ; 5 20 34 ]" << endl;
	}
	
	{
		TESTTITLE("Matrix <-> Matrix LIL conversions");
		
		SparseLIL matLIL(3, 4);
	
		matLIL(0, 1) = 2.0;
		matLIL(0, 3) = 4.0;
		matLIL(1, 1) = 1.0;
		matLIL(2, 1) = 5.0;
		matLIL(2, 0) = 3.0;
		
		DUMP_MAT(matLIL);
		
		Sparse m(matLIL);
		DUMP_MAT(m);
		
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
		DUMP_MAT(s);
		
		SparseLIL nLIL(s);
		DUMP_MAT(nLIL);
	}
	
	
	{
		TESTTITLE("Jacobi solver");
	
		// [ 2 1 ; 5 7 ] * [ x ; y ] = [ 11 ; 13 ]

		SparseLIL sLIL(2, 2);
		sLIL(0, 0) = 2.0;
		sLIL(0, 1) = 1.0;
		sLIL(1, 0) = 5.0;
		sLIL(1, 1) = 7.0;
		
		Sparse s(sLIL);
		DUMP_MAT(s);

		Vector b(2);
		b(0) = 11.;
		b(1) = 13.;

		Vector solution = s.jacobi(b);
		DUMP_VEC(solution);
		
		cout << "Expected: [ 7.1111 ; -3.2222 ]" << endl;
	}

	{
		TESTTITLE("Conjugate Gradient solver");
	
		{
			// [ 4 1 ; 1 3 ] * [ 0.0909 ; 0.6364 ] = [ 1 ; 2 ]
			
			SparseLIL sLIL(2, 2);
			sLIL(0, 0) = 4.0;
			sLIL(0, 1) = 1.0;
			sLIL(1, 0) = 1.0;
			sLIL(1, 1) = 3.0;
			
			Sparse s(sLIL);
			DUMP_MAT(s);

			Vector b(2);
			b(0) = 1.0;
			b(1) = 2.0;
			
			Vector solution = s.conjGradient(b);
			DUMP_VEC(solution);
			
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
			DUMP_MAT(s);

			Vector b(3);
			b(0) = 1.0;
			b(1) = 2.0;
			b(2) = 3.0;
			
			Vector solution = s.conjGradient(b);
			DUMP_VEC(solution);
			
			cout << "Expected: [ 2.5 ; 4 ; 3.5 ]" << endl;
		}
	}
	
	{
		TESTTITLE("LU solver");
	
		{
			// [ 4 1 ; 1 3 ] * [ 0.0909 ; 0.6364 ] = [ 1 ; 2 ]
			
			SparseLIL sLIL(2, 2);
			sLIL(0, 0) = 4.0;
			sLIL(0, 1) = 1.0;
			sLIL(1, 0) = 1.0;
			sLIL(1, 1) = 3.0;
			
			Sparse s(sLIL);
			DUMP_MAT(s);

			Vector b(2);
			b(0) = 1.0;
			b(1) = 2.0;
			
			Vector solution = s.LU(b);
			DUMP_VEC(solution);
			
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
			DUMP_MAT(s);

			Vector b(3);
			b(0) = 1.0;
			b(1) = 2.0;
			b(2) = 3.0;
			
			Vector solution = s.LU(b);
			DUMP_VEC(solution);
			
			cout << "Expected: [ 2.5 ; 4 ; 3.5 ]" << endl;
		}
		
		{
			// example where the pivot gets zero
			// taken from G. Allaire, S. M. Kaber: "Numerical Linear Algebra"
			// [ 2 4 -4 1 ; 3 6 1 -2 ; -1 1 2 3 ; 1 1 -4 1 ] * [ 1 ; -1 ; 0 ; 2 ] = [ 0 ; -7 ; 4 ; 2 ]
			
			SparseLIL sLIL(4, 4);
			sLIL(0, 0) = 2.0;
			sLIL(0, 1) = 4.0;
			sLIL(0, 2) = -4.0;
			sLIL(0, 3) = 1.0;
			sLIL(1, 0) = 3.0;
			sLIL(1, 1) = 6.0;
			sLIL(1, 2) = 1.0;
			sLIL(1, 3) = -2.0;
			sLIL(2, 0) = -1.0;
			sLIL(2, 1) = 1.0;
			sLIL(2, 2) = 2.0;
			sLIL(2, 3) = 3.0;
			sLIL(3, 0) = 1.0;
			sLIL(3, 1) = 1.0;
			sLIL(3, 2) = -4.0;
			sLIL(3, 3) = 1.0;
			
			Sparse s(sLIL);
			DUMP_MAT(s);

			Vector b(4);
			b(0) = 0.0;
			b(1) = -7.0;
			b(2) = 4.0;
			b(3) = 2.0;
			
			Vector solution = s.LU(b);
			DUMP_VEC(solution);
			
			cout << "Expected: [ 1 ; -1 ; 0 ; 2 ]" << endl;
		}
	}

	{
		TESTTITLE("Matrix load");
	
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		DUMP_MAT(m);
	}

	{
		TESTTITLE("Matrix multiplication");
	
		Sparse m;
		ifstream f("data/mtest01.dat");
		f >> m;
		
		Vector v(3);
		v(0) = 1.;
		v(1) = 2.;
		v(2) = 3.;
		
		Vector w = m * v;
		
		DUMP_VEC(w);
	}

	{
		TESTTITLE("Matrix product v' * M * w");
		
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
		TESTTITLE("Mesh load");
		
        // ifstream meshfile;
        // meshfile.open("data/mesh/mesh.msh", ios::in);
        
        Mesh M("data/mesh/mesh.msh");
    }

    {
		TESTTITLE("Mesh load");
		
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
		
		//return 0;
		
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
    	TESTTITLE("Basic triangle testing = ");

    	// construct mesh with 3 vertices, 1 triangle and one edge
    	Mesh M(3,1,3); 

    	Vertex A(0, 0, 0, 0);
    	Vertex B(0, 1, 0, 1);
    	Vertex C(1, 0, 0, 2);

    	M.V[0] = A;
    	M.V[1] = B;
    	M.V[2] = C;

    	Triangle T(&M.V[0], &M.V[1], &M.V[2], 0, 0);

    	M.T[0] = T;


    	M.V[0].T.push_back(&M.T[0]);
    	M.V[1].T.push_back(&M.T[0]);
    	M.V[2].T.push_back(&M.T[0]);

    	cout << *M.V[0].T[0];

    	BoundEdge AB(&M.V[0], &M.V[1], 0, 0);
    	BoundEdge AC(&M.V[0], &M.V[2], 0, 1);
    	BoundEdge BC(&M.V[1], &M.V[2], 0, 2);


        M.E[0] = AB;
        M.E[1] = AC;
        M.E[2] = BC;

    	Vector nAB(2);
    	nAB.constructNormal(M.E[0]);

    	Vector nAC(2);
    	nAC.constructNormal(M.E[1]);

    	Vector nBC(2);
    	nBC.constructNormal(M.E[2]);


    	cout << M.E[0];
    	cout << " Normal: " << nAB(0) << " " << nAB(1) << endl;

    	cout << M.E[1];
    	cout << " Normal: " << nAC(0) << " " << nAC(1) << endl;

    	cout << M.E[2];
    	cout << " Normal: " << nBC(0) << " " << nBC(1) << endl;
    }
	
	{
		TESTTITLE("Linear combinaton of base functions");
		
        Mesh mesh("data/mesh/square_9.msh");
		
		Vector uh(mesh.countVertices());
		
		for (uint i = 0; i < uh.size(); i++)
		{
			uh(i) = i;
		}
		
		cout << "(0.3, 0.3): " << mesh.eval(0.3, 0.3, uh);
		cout << "(0.5, 0.3): " << mesh.eval(0.5, 0.3, uh);
		cout << "(0.1, 0.7): " << mesh.eval(0.1, 0.7, uh);
		
		// points outside of the mesh
		mesh.eval(1.5, 0.3, uh);
		mesh.eval(0.3, 1.5, uh);
    }

	{
		TESTTITLE("Plot linear combination of  base functions");
		
        Mesh mesh("data/mesh/square_9.msh");
		
		Vector uh(mesh.countVertices());
		
		for (uint i = 0; i < uh.size(); i++)
		{
			uh(i) = i;
		}
		
		PlotMesh plot("test_eval", mesh, uh, "Combination of base functions");
		plot.generate(ePT_GNUPLOT_SURF, true);
	}
	
	return 0;
}
