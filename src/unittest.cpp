/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: unittest.cpp

 Description: test

**************************************************************/

#include <iostream>

#include "algebra.h"
#include "mesh.h"
#include "visualization.h"

using namespace std;
/*
// @@@ when are const/nonconst values used as rvalue/lvalue?
class _CHelper2
{
private:
	double &value;

public:
	_CHelper2(double &v) : value(v)
	{
	}
	
	double& operator= (const double& v) 
	{
		cout << "as lvalue" << endl;
		value = v;
		return value;
	}
	
	operator const double& () const
	{
		cout << "as rvalue" << endl;
		return value;
	}
};

class CTest
{
private:
	double vals[5];

public:
	_CHelper2 operator[] (uint i)
	{
		cout << "nonconst ";
		return _CHelper2(vals[i]);
	}

	//double operator[] (uint i) const
	const double& operator[] (uint i) const
	{
		cout << "const as rvalue" << endl;
		return vals[i];
	}
};*/

class CTestSimple
{
private:
	double vals[3];

public:
	double& operator[] (int i)
	{
		cout << "nonconst" << endl;
		return vals[i];
	}

	double operator[] (int i) const
	{
		cout << "const" << endl;
		return vals[i];
	}
	
	double& get(int i)
	{
		cout << "get" << endl;
		return vals[i];
	}
};


// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	cout << "== Unittest ==" << endl;
/*
	// shows problems with evaluation rValue vs lValue at const/non-const
	CTest T;
	cout << T[0] << " " << T[1] << endl;
	double d = T[2];
	double d2 = T[3] + T[4];
	
	T[0] = 1.0;
	cout << T[0] << " " << T[1] << endl;
	
	return 0;*/
	
	/*CTestSimple T;
	T[0] = 1.0;
	T[1] = 2.0;
	T[2] = 3.0;
	double d = T[2];
	cout << T[0] + T[1] << endl;
	T.get(0) = 2.0;
	cout << T[0] + T[1] << endl;
	return 0;*/
	
	// ---
	
	// Algebra
	
	{
		cout << "Default constructor" << endl;
		
		CVector v;
		
		// uncomment to get assertion fault
		//v[1] = 1;
		
		cout << v << endl;
	}

	{
		cout << "Constructor with double array, no data" << endl;
		
		CVector v(3);
		
		cout << v << endl;
		
		// ---
		
		cout << "Constructor with double array" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		
		CVector w(3, data);
		
		cout << w << endl;
	}

	{
		cout << "Copy Constructor" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		CVector v(3, data);
		
		cout << v << endl;
		
		CVector w = v;
		
		cout << w << endl;
		
		v[1] = 5.0;
		w[1] = -5.0;
		
		cout << v << endl;
		cout << w << endl;
	}
	
	{
		cout << "Size" << endl;
		
		CVector v(3);
		CVector w(5);
		CVector x;
		
		cout << v.size() << " " << w.size() << " " << x.size() << " " << endl;
	}
	
	{
		cout << "Zero" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		CVector v(3, data);
		
		cout << v << endl;
		
		v.zero();
		
		cout << v << endl;
	}
	
	{
		cout << "operator []" << endl;
		CVector v(3);
		
		// Assignment
		v[0] = -2.0;
		v[1] = 2.0;
		v[2] = 3.0;
		
		cout << v << endl;
		
		cout << v[1] << endl;
	}
	
	{
		cout << "Unary operator -" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		CVector v(3, data);
		
		cout << -v << endl;
	}

	{
		cout << "Operators += and -=" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		double data2[] = {1.0, 1.0, 1.0};
		CVector v(3, data);
		CVector w(3, data2);
		
		v += w;
		cout << v << endl;
		
		CVector x(3, data);
		x -= w;
		cout << x << endl;
		
		// assertion faults
		//CVector y(2);
		//x += y;
		//x -= y;
	}
	
	{
		cout << "Operators + and -" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		double data2[] = {1.0, 1.0, 1.0};
		CVector v(3, data);
		CVector w(3, data2);
		
		cout << v + w << endl;
		cout << v - w << endl;
		
		// assertion faults
		//CVector y(2);
		//v + y;
		//v - y;
	}
	
	{
		cout << "Scalar multiplication" << endl;
		
		double data[] = {1.0, 2.0, 3.0};
		CVector v(3, data);
		cout << 2.0 * v << endl;
		
		double data2[] = {1.0, 2.0};
		CVector w(2, data2);
		cout << w * -3.0 << endl;
		
		w *= 2.0;
		cout << w * -3.0 << endl;
	}
	
	// ---
	
	{
		cout << "Matrix Dense" << endl;
	
		double data[] = {1.0, 2.0, 3.0, 4.0};
		
		CMatrix m(2, 2, eMDT_DENSE, data);
		
		cout << m << endl;
		
		m._l(0, 1) = 5.0;
		m._l(1, 0) = -1.0;
		
		cout << m << endl;
		
		cout << "Copy" << endl;
		CMatrix n = m;
		cout << n << endl;
	}
	
	{
		cout << "Matrix Diag" << endl;
	
		double data[] = {1.1, 2.2, 3.3, 4.4};
		
		CMatrix m(4, 4, eMDT_DIAG, data);
		
		cout << m << endl;
		
		m._l(2, 2) = 12.0;
		m._l(3, 3) = -1.0;
		
		cout << m << endl;
		
		cout << "Copy" << endl;
		CMatrix n = m;
		cout << n << endl;
		
		// assertion faults
		//m._l(0, 1) = 5.0;
		//m._l(1, 0) = -1.0;
	}
	
	{
		cout << "Matrix Sparse" << endl;
	
		/*CMatrix m(4, 4, eMDT_SPARSE);
		
		cout << m << endl;
		
		m._l(2, 2) = 12.0;
		m._l(3, 3) = -1.0;
		
		cout << m << endl;
		
		cout << "Copy" << endl;
		CMatrix n = m;
		cout << n << endl;*/
	}
	
	{
		cout << "Matrix Operations" << endl;
		
		double data[] = {1.0, 2.0, 3.0, 4.0};
		
		CMatrix m(2, 2, eMDT_DENSE, data);
		
		cout << m << endl;
		
		cout << -m << endl;
		
		double data2[] = {5.0, 4.0, 3.0, 2.0};
		CMatrix n(2, 2, eMDT_DENSE, data2);
		
		double data3[] = {5.0, 4.0};
		CMatrix x(1, 2, eMDT_DENSE, data3);
		
		m += n;
		cout << m << endl;
		
		// assertion fault
		//m += x;
		
		m -= n;
		cout << m << endl;
		
		// assertion fault
		//m -= x;
		
		cout << m + n << endl;
		
		// assertion fault
		//m + x;
		
		cout << m - n << endl;
		
		// assertion fault
		//m - x;
		
		double s = 1.5;
		m *= s;
		cout << m << endl;
		
		double dataVec[] = {1.0, 2.0};
		CVector v(2, dataVec);
		
		/*cout << v * m << endl;
		
		cout << m * v << endl;
		
		cout << m * s << endl;
		
		cout << s * m << endl;*/
	}
	
	// ---
	
	// Mesh
	
	{
		CMesh msh;
	}
	
	// ---
	
	// Visualization
	
	{
		CPlot plot;
	}
	
	return 0;
}