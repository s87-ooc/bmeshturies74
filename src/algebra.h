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

using namespace std;

typedef unsigned int uint;

/** values with module below are treated as zero */
const double ZERO_EPSILON = 0.00000001;

// Note that depending on the constructor call the values aren't initialized

class CVector
{
private:
	double* m_Entries;
	uint m_Size;
	
public:
	CVector();
	~CVector();
	
	/** creates a vector with "size" entries, filled with data if it's not 0 */
	CVector(uint size, double* data = 0);
	
	/** copy constructor */
	CVector(const CVector& v);
	
	/** returns the number of entries */
	uint size() const;
	
	/** set all components to 0 */
	void zero();

	/** get the i-th component */
	const double& operator[] (uint i) const;
	
	/** access the i-th component for writing */
	double& operator[] (uint i);
	
	/** unary operators */
	CVector operator- ();
	
	/** vector addition, etc */
	CVector& operator+= (const CVector& v);
	CVector& operator-= (const CVector& v);
	CVector operator+ (const CVector& v);
	CVector operator- (const CVector& v);
	
	/** scalar multiplication */
	CVector& operator*= (const double s);
};

/** scalar multiplications */
CVector operator* (const double s, const CVector& v);
CVector operator* (const CVector& v, const double s);

// operations with matrices defined below

ostream& operator<<(ostream& o, const CVector& v);
//istream& operator>>(istream& i, CVector& v);

// ----------------------------------------------------------------------------

// how data is stored in a matrix depends on what type of matrix we have

class CMatrixData
{
friend class CMatrix;
protected:
	uint m_SizeRows;
	uint m_SizeColumns;

public:
	static const double s_Zero;

	/** access the (r, c) component non-const (rValue/lValue) */
	//virtual double& operator() (uint r, uint c) = 0;

	/** access the (r, c) component const (rValue only) */
	//virtual const double& operator() (uint r, uint c) const = 0;
	
	// @ workaround
	virtual double& _l(uint r, uint c) = 0;
	virtual const double& _r(uint r, uint c) const = 0;
	// @@@
};

// ----------------------------------------------------------------------------

class CMatrixDataDense : public CMatrixData
{
private:
	double* m_Data;

public:
	CMatrixDataDense();
	~CMatrixDataDense();

	CMatrixDataDense(uint rows, uint columns, double* data = 0);

	/** copy constructor */
	CMatrixDataDense(const CMatrixDataDense& mdata);
	
	//double& operator() (uint r, uint c);
	//const double& operator() (uint r, uint c) const;
	
	// @ workaround
	double& _l(uint r, uint c);
	const double& _r(uint r, uint c) const;
	// @@@
};

// ----------------------------------------------------------------------------

class CMatrixDataDiag : public CMatrixData
{
private:
	double* m_Data;

public:
	CMatrixDataDiag();
	~CMatrixDataDiag();

	CMatrixDataDiag(uint size, double* data = 0);

	/** copy constructor */
	CMatrixDataDiag(const CMatrixDataDiag& mdata);
	
	//double& operator() (uint r, uint c);
	//const double& operator() (uint r, uint c) const;
	
	// @ workaround
	double& _l(uint r, uint c);
	const double& _r(uint r, uint c) const;
	// @@@
};

// ----------------------------------------------------------------------------

class CMatrixDataSparse : public CMatrixData
{
private:
	double* m_Data;

public:
	CMatrixDataSparse();
	~CMatrixDataSparse();

	CMatrixDataSparse(uint rows, uint columns);
	
	/** copy constructor */
	CMatrixDataSparse(const CMatrixDataSparse& mdata);
	
	//double& operator() (uint r, uint c);
	//const double& operator() (uint r, uint c) const;
	
	// @ workaround
	double& _l(uint r, uint c);
	const double& _r(uint r, uint c) const;
	// @@@
};

// ----------------------------------------------------------------------------

enum eMatrixDataType
{
	eMDT_DENSE = 0,
	eMDT_DIAG,
	eMDT_SPARSE
};

// Note that depending on the matrix data type the entries aren't initialized

class CMatrix
{
private:
	/** data depends on the type of the matrix */
	CMatrixData* m_Data;
	
	/** stores the type of the matrix */
	eMatrixDataType m_Type;
	
public:
	CMatrix();
	~CMatrix();
	
	/** general constructor for matrices, data is only used for DENSE and DIAG matrices  */
	CMatrix(uint rows, uint columns, eMatrixDataType type = eMDT_DENSE, double* data = 0);
	
	/** copy constructor */
	CMatrix(const CMatrix& m);
	
	/** returns the number of rows */
	uint sizeRows() const;
	
	/** returns the number of columns */
	uint sizeColumns() const;

	//double& operator() (uint r, uint c);	
	//const double& operator() (uint r, uint c) const;
	
	// @ workaround
	double& _l(uint r, uint c);
	const double& _r(uint r, uint c) const;
	// @@@
	
	/** unary operators */
	CMatrix operator- ();
	
	/** matrix addition, etc */
	CMatrix& operator+= (const CMatrix& m);
	CMatrix& operator-= (const CMatrix& m);
	CMatrix operator+ (const CMatrix& m);
	CMatrix operator- (const CMatrix& m);
	
	/** multiplication with scalar */
	CMatrix& operator *= (const double s);
};

/** multiplication with vector from left gives vector */
CVector operator* (const CVector& v, const CMatrix& m);
/** multiplication with vector from right gives vector */
CVector operator* (const CMatrix& m, const CVector& v);

/** multiplication with scalar from left */
CMatrix operator* (const double s, const CMatrix& m);
/** multiplication with scalar from right */
CMatrix operator* (const CMatrix& m, const double s);

ostream& operator<<(ostream& o, const CMatrix& m);

#endif // __ALGEBRA_H__