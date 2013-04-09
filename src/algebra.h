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
#include <map>

#include "types.h"


// ----------------------------------------------------------------------------

typedef std::vector<double> TVals;
typedef std::map<int, double> TValsMap;
typedef std::vector<int> TInd;

// ----------------------------------------------------------------------------

class Sparse;

class Vector
{
private:
	TVals mVals;

	// ---

	friend std::istream& operator>> (std::istream &in, Vector& v);
	friend std::ostream& operator<< (std::ostream &out, Vector& v);
	
public:
	Vector();
    Vector(uint size);
	Vector(TVals &vals);
    
	uint size() const;

	/** Euclidean norm */
	double norm2() const;

	double operator() (uint pos) const;
	double& operator() (uint pos);
	
	Vector operator+ (Vector const& v) const;
	Vector operator- (Vector const& v) const;
	
	Vector& operator+= (Vector& v);
	Vector& operator-= (Vector& v);

	Vector& operator*= (const double s);
	
	//void constructF(const Mesh& m);
	//void constructG(const Mesh& m);
};

Vector operator* (const Vector& v, const double s);
Vector operator* (const double s, const Vector& v);

/// dot product v1 * v2
double dot(const Vector& v1, const Vector& v2);

// ----------------------------------------------------------------------------

/** Sparse matrix stored in LIL (List of Lists) format
	Optimized for manipulation of single entries/rows.
	Not effective for Matrix/Vector operations, convert to Sparse instead.
	
	We're using this class for the LU solver, Sparse should be preferred.

	General assumption: we don't have empty rows */

class Sparse;

class SparseLIL
{
private:
	/** values per row */
	TVals mRowVals[]; // should we use pointers here?
	/** column indices per row */
	TInd mColInd[];
	
	uint mSizeRows;
	uint mSizeColumns;
	
	// ---
	
	friend class Sparse;

public:
	SparseLIL();
	SparseLIL(const SparseLIL& m);
	
	/** create empty sparse matrix with specified dimensions */
	SparseLIL(uint rows, uint columns);
	
	/** convert CSR to LIL */
	SparseLIL(const Sparse& matCSR);

    double operator() (uint row, uint col) const;
    double& operator() (uint row, uint col);

	uint sizeColumns() const;
    uint sizeRows() const;
    uint sizeNNZ() const;
};

// ----------------------------------------------------------------------------

/** Sparse matrix stored in CSR (Compressed Sparse Row) format
	Optimized for Matrix/Vector operations.
	After the matrix was constructed, it's no longer possible to manipulate
	elements directly, convert to SparseLIL instead.
	
	General assumption: we don't have empty rows */

class Mesh;

class Sparse
{
private:
    TVals mVals;
    TInd mColInd;
	/** number of rows + 1, last entry is last nnz in row + 1 */
	TInd mRowPtr;

	uint mSizeColumns;

	// ---

	friend class SparseLIL;
	
	/** format: "i j value", row i ascending */
    friend std::istream& operator >>(std::istream&, Sparse&);
    friend std::ostream& operator <<(std::ostream&, Sparse&);

	/// (v' * m)' ?!
	friend Vector operator* (const Vector& v, const Sparse& m);
	friend Vector operator* (const Sparse& m, const Vector& v);

	/// v1' * m * v2
	friend double prod (const Vector& v1, const Sparse& m, const Vector& v2);
    
public:
    Sparse();
    Sparse(const Sparse& m);

    Sparse(const TVals& vals, const TInd& colInd, const TInd& rowPtr, uint columns);
    Sparse(const TValsMap& M, uint rows, uint cols);
	
	/** convert LIL to CSR */
	Sparse(const SparseLIL& matLIL);

    uint sizeColumns() const;
    uint sizeRows() const;
    uint sizeNNZ() const;

	/** values can only be accessed read-only  */
    double operator() (uint row, uint col) const;
    //double& operator() (uint row, uint col);
    
	/** solve Ax = b with jacobi */
    Vector jacobi(Vector const& v) const;

	/** solve Ax = b with conjugate Gradient */
    Vector conjGradient(Vector const& v) const;

	/** solve Ax = b with LU, return decomposition if needed */
    Vector LU(Vector const& v, Sparse* m = 0) const;
	
	void constructA(const Mesh& m);
	void constructM(const Mesh& m);
	void constructB(const Mesh& m);
};

#endif // __ALGEBRA_H__
