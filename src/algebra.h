#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenković (stjepan@stjepan.net)

 ---

     Fichier: algebra.h

 Description: déclaration des classes mathématiques pour les matrices, vecteurs, etc

**************************************************************/

#include <iostream>

#include <vector>
#include <map>

#include "types.h"
#include "mesh.h"

// ----------------------------------------------------------------------------

#define JACOBI_TOLERANCE 0.000000001
#define CONJGRADIENT_TOLERANCE 0.000000001

// ----------------------------------------------------------------------------

typedef std::vector<double> TVals;
typedef std::map<int, double> TValsMap;
typedef std::vector<uint> TInd;

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

	double operator() (Vertex* v) const;
	
	Vector operator+ (Vector const& v) const;
	Vector operator- (Vector const& v) const;
	
	Vector& operator+= (Vector& v);
	Vector& operator-= (Vector& v);

	Vector& operator*= (const double s);
	
	/** construct the vector based on vertex values evaluated in f */
	Vector& constructFunc(const Mesh& mesh, double (*f)(const Vertex&));
	
	/** construct the vector based on integration over the edges of the mesh */
	Vector& constructFuncIntSurf(const Mesh& mesh, double (*f)(const Vertex&));

	/** construct the vector based on integration over the edges of the mesh 
	taking into account the normal vector to the domain boundary */
	Vector& constructFuncSurf(const Mesh& mesh, double (*f)(const Vertex&, const BoundEdge&));

	/** construct normal vector to edge */
	Vector& constructNormal(const BoundEdge& edge);

};

Vector operator* (const Vector& v, const double s);
Vector operator* (const double s, const Vector& v);

/// dot product v1 * v2
double dot(const Vector& v1, const Vector& v2);

double globalL2Error(const Mesh& mesh, const Vector& exact, const Vector& approx);


// ----------------------------------------------------------------------------

/** Sparse matrix stored in LIL (List of Lists) format
	Optimized for manipulation of single entries/rows.
	Not effective for Matrix/Vector operations, convert to Sparse instead.
	
	We're using this class for the LU solver, Sparse should be preferred.

	General assumption: we don't have empty rows */

class Sparse;

class SparseMap
{
private:
	TValsMap mVals;
	uint mSizeRows;
	uint mSizeColumns;

	friend class Sparse;
	friend std::ostream& operator <<(std::ostream&, SparseMap& m);


public:
	SparseMap(uint rows, uint columns);

	double operator() (uint row, uint col) const;
	// double& operator() (uint row, uint col);
	void addAt(uint row, uint col, double value);
	void setZero(uint row, uint col);

	SparseMap& operator+=(const SparseMap& rhs);
	SparseMap& operator*=(const double& val);

	uint sizeColumns() const;
    uint sizeRows() const;
    uint sizeNNZ() const;

    /** Assumption: mVals is empty when construct is called */
	SparseMap& constructA(const Mesh& mesh);
	SparseMap& constructM(const Mesh& mesh);
	SparseMap& constructB(const Mesh& mesh);
};

class SparseLIL
{
private:
	TVals* mRowVals; /** values per row */
	TInd* mColInd; /** column indices per row */
	
	uint mSizeRows;
	uint mSizeColumns;
	
	// ---
	
	friend class Sparse;

public:
	SparseLIL();
	~SparseLIL();
	
	SparseLIL(const SparseLIL& m);
	
	/** create empty sparse matrix with specified dimensions */
	SparseLIL(uint rows, uint columns);
	
	/** convert CSR to LIL */
	SparseLIL(const Sparse& matCSR);

	/** return the product of the matrix and its transpose: M * M' */
	SparseLIL prodTranspose() const;
	
	/** get value read-only */
    double operator() (uint row, uint col) const;
	
	/** before assigning values at (i, j) make sure not to add zeros */
    double& operator() (uint row, uint col);

	uint sizeColumns() const;
    uint sizeRows() const;
    uint sizeNNZ() const;
};

// ----------------------------------------------------------------------------

/** Sparse matrix stored in CSR (Compressed Sparse Row) format
	Optimized for Matrix/Vector operations, features solvers.
	After the matrix was constructed, it's no longer possible to manipulate
	elements directly, convert to SparseLIL instead.
	
	General assumption: we don't have empty rows */

class Mesh;

class Sparse
{
private:
    TVals mVals;
    TInd mColInd;
	TInd mRowPtr; /** number of rows + 1, last entry is last nnz in row + 1 */

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

    /** convert SparseMap to CSR */
    Sparse(const SparseMap& M);
	
	/** convert LIL to CSR */
	Sparse(const SparseLIL& matLIL);

    uint sizeColumns() const;
    uint sizeRows() const;
    uint sizeNNZ() const;

	/** values can only be read  */
    double operator() (uint row, uint col) const;
	
	// ---

	/** solve Ax = b with conjugate Gradient
		Assumption: the matrix is symmetric positive-definite */
    Vector conjGradient(Vector const& b) const;
    
	/** solve Ax = b with jacobi
		Assumption: the matrix is strictly diagonally dominant AND 2D - A is too */
    Vector jacobi(Vector const& b) const;
	/** solve Ax = b with jacobi and return whether the method converged */
	Vector jacobi(Vector const& b, bool& convergence) const;

	/** solve Ax = b with LU, optionally return the decomposed matrix */
    Vector LU(Vector const& b, Sparse* lu = 0) const;
	
	/** calculate uNew and vNew using the Newmark method using this matrix as mass matrix M
		according to the notation from the course we have:
		M * uNew = uMatU * u + vMatU * v + uFactor * (F(t) + G(t))
		M * vNew = uMatV * (u + uNew) + M * v + vFactor * (F(t) + G(t) + F(t + dt) + G(t + dt)) */
	void newmark(Vector& uNew, Vector& vNew, const Vector& u, const Vector& v, const Sparse& uMatU, const Sparse& vMatU, const Sparse& uMatV, double uFactor = 0., double vFactor = 0., double dt = 0., double t = 0., double (*f)(const Vertex&, double) = 0, double (*g)(const Vertex&, double) = 0) const;
};

#endif // __ALGEBRA_H__
