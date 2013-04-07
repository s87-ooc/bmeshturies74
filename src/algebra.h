#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: algebra.h

 Description: déclaration des classes mathématiques pour les matrices, vecteurs, etc

**************************************************************/


#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

typedef vector<double> TVals;
typedef vector<int> TInd;

// ----------------------------------------------------------------------------

class Sparse;

class Vector
{
private:
	TVals mVals;

	friend istream& operator>> (istream &in, Vector& v);
	friend ostream& operator<< (ostream &out, Vector& v);
	
public:
	Vector();
    Vector(unsigned int size);
	Vector(TVals &vals);
    
	unsigned int size() const;

	/** Euclidean norm */
	double norm2() const;

	double operator() (unsigned int pos) const;
	double& operator() (unsigned int pos);
	
	Vector operator+ (Vector const& v) const;
	Vector operator- (Vector const& v) const;
	
	Vector& operator+= (Vector& v);
	Vector& operator-= (Vector& v);

	Vector& operator*= (const double s);
};

Vector operator* (const Vector& v, const double s);
Vector operator* (const double s, const Vector& v);

/// dot product v1 * v2
double dot(const Vector& v1, const Vector& v2);

// ----------------------------------------------------------------------------

class Sparse
{
private:
    TVals mVals;
    TInd mColInd, mRowPtr;
    
    friend istream& operator >>(istream&, Sparse&);
    friend ostream& operator <<(ostream&, Sparse&);

	friend Vector operator* (const Vector& v, const Sparse& m);
	friend Vector operator* (const Sparse& m, const Vector& v);

	/// x' * m * x
	friend double prod (const Vector& v, const Sparse& m);
    
public:
    Sparse();
    Sparse(const Sparse& m);

    Sparse(const TVals& vals, const TInd& colInd, const TInd& rowPtr);
    Sparse(unsigned int nrows, unsigned int ncols);
    
    unsigned int sizeColumns() const;
    unsigned int sizeRows() const;
    unsigned int sizeNNZ() const;
    
    double operator() (unsigned int col, unsigned int row) const;
    double& operator() (unsigned int col, unsigned int row);
    
    Vector jacobi(Vector const& v) const;
    Vector gradient_conj(Vector const& v) const;
};


#endif // __ALGEBRA_H__
