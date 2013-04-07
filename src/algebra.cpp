/*************************************************************

 Mini Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: algebra.cpp

 Description: implementation des classes mathématiques pour les matrices, vecteurs, etc

**************************************************************/

#include "algebra.h"

#include <cstring>
#include <assert.h>
#include <math.h>

// ----------------------------------------------------------------------------

#define JACOBI_TOLERANCE 0.000001

// ----------------------------------------------------------------------------

Vector::Vector()
{
}

Vector::Vector(unsigned int size) :
mVals(size)
{
}

Vector::Vector(TVals &vals):
mVals(vals)
{
}

unsigned int Vector::size() const
{
	return mVals.size();
}


double Vector::norm2() const
{
	return sqrt(dot(*this, *this));
}

double Vector::operator() (unsigned int pos) const
{
	return mVals[pos];
}

double& Vector::operator() (unsigned int pos)
{
	return mVals[pos];
}

Vector Vector::operator+ (Vector const& v) const
{
	assert(v.size() == mVals.size());

    int n = v.size();
	Vector res(n);
	
	for (int i = 0; i < n; i++)
		res(i) = mVals[i] + v(i);
	
	return res;
}


Vector Vector::operator- (Vector const& v) const
{
	assert(v.size() == mVals.size());

    int n = v.size();
	Vector res(n);
	
	for (int i = 0; i < n; i++)
		res(i) = mVals[i] - v(i);
	
	return res;
}


Vector& Vector::operator+= (Vector &v)
{
	assert(v.size() == mVals.size());

	for (int i = 0; i < v.size(); i++)
		mVals[i] += v(i);
    
	return *this;
}

Vector& Vector::operator-= (Vector &v)
{
	assert(v.size() == mVals.size());

	for (int i = 0; i < v.size(); i++)
		mVals[i] -= v(i);
	
	return *this;
}

Vector& Vector::operator*= (const double s)
{
	for (int i = 0; i < mVals.size(); i++)
		mVals[i] *= s;
	
	return *this;
}

istream& operator>> (istream &in, Vector& v)
{
	// TODO
	return in;
}

ostream& operator<<(ostream &out, Vector& v)
{
	cout << "[ ";
	for (int i = 0; i < v.size(); i++)
	{
		out << v(i) << " ";
	}
	cout << "] " << v.size();
	
	return out;
}

Vector operator* (const Vector& v, const double s)
{
	Vector res(v.size());

	for (int i = 0; i < res.size(); i++)
		res(i) = v(i) * s;
	
	return res;
}

Vector operator* (const double s, const Vector& v)
{
	Vector res(v.size());

	for (int i = 0; i < res.size(); i++)
		res(i) = v(i) * s;
	
	return res;
}

double dot(const Vector& v1, const Vector& v2)
{
	assert(v1.size() == v2.size());

	double res = 0.0;

	for (int i = 0; i < v1.size(); i++)
		res += v1(i) * v2(i);

	return res;
}

// ----------------------------------------------------------------------------

Sparse::Sparse()
{
}

Sparse::Sparse(const Sparse& m) :
mVals(m.mVals),
mColInd(m.mColInd),
mRowPtr(m.mRowPtr),
mColumns(m.mColumns)
{
}

Sparse::Sparse(const TVals& vals, const TInd& colInd, const TInd& rowPtr, unsigned int columns) :
mVals(vals),
mColInd(colInd),
mRowPtr(rowPtr),
mColumns(columns)
{
}

unsigned int Sparse::sizeColumns() const
{
    return mColumns;
}

unsigned int Sparse::sizeRows() const
{
	assert(mRowPtr.size() > 0);

    return mRowPtr.size() - 1;
}

unsigned int Sparse::sizeNNZ() const

{
    return mVals.size();
}

double Sparse::operator() (unsigned int row, unsigned int col) const
{
    assert(row < sizeRows() && col < sizeColumns());

    for (unsigned int i = mRowPtr[row]; i < mRowPtr[row+1]; i++)
    {
        if (mColInd[i] == col)
			return mVals[i];
    }
    
    return 0.0;
}

double& Sparse::operator() (unsigned int row, unsigned int col)
{
    assert(row < sizeRows() && col < sizeColumns());

    for (unsigned int i = mRowPtr[row]; i < mRowPtr[row+1]; i++)
    {
        if (mColInd[i] == col)
			return mVals[i];
    }
    
	// we never add new nnz elements at this point
	cout << "Can't add new NNZ element here!" << endl;
	assert(false);

	// won't be called
    return mVals[0];
}

Vector Sparse::jacobi(Vector const& b) const
{
    Vector x1(sizeColumns());
    Vector x2(sizeColumns());
    Vector& xcur = x1;
    Vector& xnew = x2;
    
    unsigned int k = 0;
    double diag;
    bool convergence = false;
    
    while (!convergence)
    {
        for(unsigned int i = 0; i < sizeRows(); i++)
        {
            xnew(i) = b(i);
            diag = 0;
            // iterate over the ith row of the matrix

            for (k = mRowPtr[i]; k < mRowPtr[i+1]; k++)
            {
                if(mColInd[k] != i)
					xnew(i) -= mVals[k] * xcur(mColInd[k]);

                else
                    diag = mVals[k];
            }

            xnew(i) /= diag;
        }

		convergence = (xnew - xcur).norm2() < JACOBI_TOLERANCE;

        // implement equality operator for Vector
		Vector& xtmp = xcur;
        xcur = xnew;
		xnew = xtmp;
    }
   
    return xnew;
}

Vector Sparse::conjGradient(Vector const& b) const
{
	// TODO: Implement!
	Vector p;
    return p;
}

Vector Sparse::LU(Vector const& v, Sparse* m) const
{
	Vector p;
	return p;
}

istream& operator >> (istream &is, Sparse &m)
{
	unsigned int sizeRows, sizeColumns, sizeNNZ;
	is >> sizeRows >> sizeColumns >> sizeNNZ;

    TVals vals(sizeNNZ);    
	TInd colInd(sizeNNZ);
	TInd rowPtr(sizeRows + 1);
    
    unsigned int r;
    
    for(unsigned int k = 0; k < sizeNNZ; k++)
    {
        is >> r >> colInd[k] >> vals[k];
        if (k == 0 || r > rowPtr[k-1])
        {
			rowPtr[r] = k;
		}
    }
    
    rowPtr[sizeRows] = rowPtr[sizeRows - 1] + 1;
    
    m = Sparse(vals, colInd, rowPtr, sizeColumns);
    
    return is;
}

ostream& operator<< (ostream &os, Sparse &m)
{
	cout << m.sizeRows() << " "
		<< m.sizeColumns() << " "
		<< m.sizeNNZ() << endl;
	
	for (unsigned int i = 0; i < m.sizeRows(); i++)
	{
		for (unsigned int k = m.mRowPtr[i]; k < m.mRowPtr[i+1]; k++)
		{
			cout << i << " " << m.mColInd[k] << " " << m.mVals[k] << endl;
		}
	}

    return os;
}

Vector operator* (const Vector& v, const Sparse& m)
{
	Vector res(m.sizeColumns());
    
    // TODO (maybe)
    
    return res;
}

Vector operator* (const Sparse& m, const Vector& v)
{
    assert(m.sizeColumns() == v.size());
    
    Vector res(m.sizeColumns());
    
    for(unsigned int i = 0; i < m.sizeRows(); i++)
    {
        for(unsigned int k = m.mRowPtr[i]; k < m.mRowPtr[i+1]; k++)
        {
            res(i) += m.mVals[k] * v(m.mColInd[k]);
        }
    }
    
    return res;
}

double prod (const Vector& v1, const Sparse& m, const Vector& v2)
{
	double res;
	
	Vector vTmp = m * v2;
	res = dot(v1, vTmp);
	
	return res;
}

