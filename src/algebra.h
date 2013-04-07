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

using namespace std;

typedef vector<double> TVals;
typedef vector<int> TInd;

// ----------------------------------------------------------------------------

class Sparse;

class Vector
{
private:
	TVals mVals;

	// ---

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

/** General assumption: we don't have empty rows */

class Sparse
{
private:
    TVals mVals;
    TInd mColInd;
	/** number of rows + 1, last entry is last nnz in row + 1 */
	TInd mRowPtr;

	unsigned int mColumns;

	// ---

	/** format: "i j value", row i ascending */
    friend istream& operator >>(istream&, Sparse&);
    friend ostream& operator <<(ostream&, Sparse&);

	/// (v' * m)' ?!
	friend Vector operator* (const Vector& v, const Sparse& m);
	friend Vector operator* (const Sparse& m, const Vector& v);

	/// v1' * m * v2
	friend double prod (const Vector& v1, const Sparse& m, const Vector& v2);
    
public:
    Sparse();
    Sparse(const Sparse& m);

    Sparse(const TVals& vals, const TInd& colInd, const TInd& rowPtr, unsigned int columns);
    
    unsigned int sizeColumns() const;
    unsigned int sizeRows() const;
    unsigned int sizeNNZ() const;

    double operator() (unsigned int row, unsigned int col) const;
    double& operator() (unsigned int row, unsigned int col);
    
	/** solve Ax = b with jacobi */
    Vector jacobi(Vector const& v) const;

	/** solve Ax = b with conjugate Gradient */
    Vector conjGradient(Vector const& v) const;

	/** solve Ax = b with LU, return decomposition if needed */
    Vector LU(Vector const& v, Sparse* m = 0) const;
};


#endif // __ALGEBRA_H__
