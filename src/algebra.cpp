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

using namespace std;

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
	
	for (uint i = 0; i < n; i++)
		res(i) = mVals[i] - v(i);
	
	return res;
}


Vector& Vector::operator+= (Vector &v)
{
	assert(v.size() == mVals.size());

	for (uint i = 0; i < v.size(); i++)
		mVals[i] += v(i);
    
	return *this;
}

Vector& Vector::operator-= (Vector &v)
{
	assert(v.size() == mVals.size());

	for (uint i = 0; i < v.size(); i++)
		mVals[i] -= v(i);
	
	return *this;
}

Vector& Vector::operator*= (const double s)
{
	for (uint i = 0; i < mVals.size(); i++)
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
	for (uint i = 0; i < v.size(); i++)
	{
		out << v(i) << " ";
	}
	cout << "] " << v.size();
	
	return out;
}

Vector operator* (const Vector& v, const double s)
{
	Vector res(v.size());

	for (uint i = 0; i < res.size(); i++)
		res(i) = v(i) * s;
	
	return res;
}

Vector operator* (const double s, const Vector& v)
{
	Vector res(v.size());

	for (uint i = 0; i < res.size(); i++)
		res(i) = v(i) * s;
	
	return res;
}

double dot(const Vector& v1, const Vector& v2)
{
	assert(v1.size() == v2.size());

	double res = 0.0;

	for (uint i = 0; i < v1.size(); i++)
		res += v1(i) * v2(i);

	return res;
}

// ----------------------------------------------------------------------------

SparseMap::SparseMap(uint rows, uint columns) :
mSizeRows(rows),
mSizeColumns(columns),
mVals()
{
}

double SparseMap::operator() (uint row, uint col) const
{
	// check if there is a non-zero value at a given index
	try
	{
		return mVals.at(row * mSizeColumns + col);
	}
	catch(...)
	{
		return 0.0;
	}
}

double& SparseMap::operator() (uint row, uint col)
{
	return mVals[row*mSizeColumns + col];
}

void SparseMap::addAt(uint row, uint col, double value)
{
	// only add non-zero values to avoid filling with zeros
	if (fabs(value) > EQ_TOL)
	{
		mVals[row*mSizeColumns + col] += value;
	}
}

void SparseMap::setZero(uint row, uint col)
{
	mVals.erase(row*mSizeColumns + col);
}

SparseMap& SparseMap::operator+=(const SparseMap& rhs)
{
	// Dimension check
	assert( rhs.sizeColumns() == this->sizeColumns());
	assert( rhs.sizeRows() == this->sizeRows());

	for(TValsMap::const_iterator iter = rhs.mVals.begin(); iter != rhs.mVals.end(); iter++)
	{
		this->mVals[iter->first] += iter->second;

		if( fabs(this->mVals[iter->first]) < EQ_TOL)
		{
			this->mVals.erase(iter->first);
		}
	}
	return (*this);
}
SparseMap& SparseMap::operator*=(const double& val)
{
	for(TValsMap::iterator iter = this->mVals.begin(); iter != this->mVals.end(); iter++)
	{
		iter->second *= val;

		if( fabs(iter->second) < EQ_TOL)
		{
			mVals.erase(iter->first);
		}
	}
	return (*this);
}

uint SparseMap::sizeColumns() const
{
	return mSizeColumns;
}

uint SparseMap::sizeRows() const
{
	return mSizeRows;
}

uint SparseMap::sizeNNZ() const
{
	return mVals.size();
}

SparseMap constructM(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	SparseMap M(Nv, Nv);

	double area;
	double mel = 1/12.0;
	double value;
	int coef;

	// iterate over the triangles of the mesh
	for(  uint t = 0; t < mesh.countTriangles(); t++)
	{
		area = mesh.T[t].area;

		// iterate over vertices of the current triangle
		for(uint i=0; i < 2; i++)
		{
			for(uint j=0; j<2; j++)
			{
				coef = ( i == j ) ? 2 : 1;
				value = coef * mel * area;
				M.addAt(i, j, value);
			}
		}
	}

	return M;

}
SparseMap constructM(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	SparseMap M(Nv, Nv);

	return M;
}

SparseMap constructB(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	SparseMap M(Nv, Nv);

	return M;
}


// ----------------------------------------------------------------------------

SparseLIL::SparseLIL() :
mSizeRows(0),
mSizeColumns(0),
mRowVals(0),
mColInd(0)
{
}

SparseLIL::~SparseLIL()
{
	if (mRowVals)
	{
		delete[] mRowVals;
	}
	
	if (mColInd)
	{
		delete[] mColInd;
	}
}

SparseLIL::SparseLIL(const SparseLIL& m) :
mSizeRows(m.sizeRows()),
mSizeColumns(m.sizeColumns())
{
	mRowVals = new TVals[mSizeRows];
	mColInd = new TInd[mSizeRows];
	
	for (uint i = 0; i < mSizeRows; i++)
	{
		mRowVals[i] = m.mRowVals[i];		
		mColInd[i] = m.mColInd[i];
	}
}

SparseLIL::SparseLIL(uint rows, uint columns) :
mSizeRows(rows),
mSizeColumns(columns)
{
	mRowVals = new TVals[rows];
	mColInd = new TInd[rows];

	// all values are zero by default
}

SparseLIL::SparseLIL(const Sparse& matCSR) :
mSizeRows(matCSR.sizeRows()),
mSizeColumns(matCSR.sizeColumns())
{
	mRowVals = new TVals[mSizeRows];
	mColInd = new TInd[mSizeRows];

	for(uint i = 0; i < matCSR.sizeRows(); i++)
	{
		// K = index of values
		for (uint k = matCSR.mRowPtr[i]; k < matCSR.mRowPtr[i+1]; k++)
		{
			mRowVals[i].push_back(matCSR.mVals[k]);
			mColInd[i].push_back(matCSR.mColInd[k]);
		}
	}
}

double SparseLIL::operator() (uint row, uint col) const
{
	for (uint i = 0; i < mRowVals[row].size(); i++)
	{
		if (mColInd[row][i] == col)
		{
			return mRowVals[row][i];
		}
		else if (mColInd[row][i] > col)
		{
			return 0.0;
		}
	}
	
	return 0.0;
}

double& SparseLIL::operator() (uint row, uint col)
{
	TVals::iterator it;
	TInd::iterator itInd;
	
	for (it = mRowVals[row].begin(), itInd = mColInd[row].begin();
		it != mRowVals[row].end(), itInd != mColInd[row].end();
		++it, ++itInd)
	{
		if (*itInd == col)
		{
			return *it;
		}
		else if (*itInd > col)
		{
			break;
		}
	}
	
	// create a new entry
	TVals::iterator itNew = mRowVals[row].insert(it, 0.0);
	mColInd[row].insert(itInd, col);

	return *itNew;
}

uint SparseLIL::sizeColumns() const
{
	return mSizeColumns;
}

uint SparseLIL::sizeRows() const
{
	return mSizeRows;
}

uint SparseLIL::sizeNNZ() const
{
	uint nnz = 0;
	
	for (uint i = 0; i < mSizeRows; i++)
	{
		nnz += mRowVals[i].size();
	}
	
	return nnz;
}

// ----------------------------------------------------------------------------

Sparse::Sparse() :
mSizeColumns(0)
{
}

Sparse::Sparse(const Sparse& m) :
mVals(m.mVals),
mColInd(m.mColInd),
mRowPtr(m.mRowPtr),
mSizeColumns(m.mSizeColumns)
{
}

Sparse::Sparse(const TVals& vals, const TInd& colInd, const TInd& rowPtr, unsigned int columns) :
mVals(vals),
mColInd(colInd),
mRowPtr(rowPtr),
mSizeColumns(columns)
{
}

Sparse::Sparse(const SparseMap& M):
mVals(M.sizeNNZ()), 
mColInd(M.sizeNNZ()), 
mRowPtr(M.sizeRows()+1), 
mSizeColumns(M.sizeColumns())
{
	int k = 0;
	int index, row_cur, row_prev=-1;
	for(TValsMap::const_iterator iter = M.mVals.begin(); iter != M.mVals.end(); ++iter)
	{
		index = iter->first;
		mColInd[k] = index % mSizeColumns;
		mVals[k] = iter->second;
		row_cur = index / mSizeColumns;

		if (row_cur > row_prev)
        {
			mRowPtr[row_cur] = k;
		}
		row_prev = row_cur;
		k++;
	}

	mRowPtr[M.sizeRows()] = mRowPtr[M.sizeRows() - 1] + 1;
}

Sparse::Sparse(const SparseLIL& matLIL) :
mSizeColumns(matLIL.sizeColumns())
{
	uint k = 0;
	
	mRowPtr.push_back(0);
	
	for (uint i = 0; i < matLIL.sizeRows(); i++)
	{
		for (uint j = 0; j < matLIL.mRowVals[i].size(); j++)
		{
			mVals.push_back(matLIL.mRowVals[i][j]);
			mColInd.push_back(matLIL.mColInd[i][j]);
			k++;
		}
		mRowPtr.push_back(k);
	}
}

unsigned int Sparse::sizeColumns() const
{
    return mSizeColumns;
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
        for(uint i = 0; i < sizeRows(); i++)
        {
            // iterate over the ith row of the matrix

			xnew(i) = b(i);
            diag = 0;

            for (k = mRowPtr[i]; k < mRowPtr[i+1]; k++)
            {
                if(mColInd[k] != i)
					xnew(i) -= mVals[k] * xcur(mColInd[k]);

                else
                    diag = mVals[k];
            }

            xnew(i) /= diag;
        }

		// Stjepan: the convergence criterium may not be correct in this way,
		// we should do either |Axnew - b| < epsilon (costly) or
		// using the "vecteur residu" as described on http://fr.wikipedia.org/wiki/M%C3%A9thode_de_Jacobi
		// TODO: discuss convergence criterium
		convergence = (xnew - xcur).norm2() < JACOBI_TOLERANCE;

		Vector& xtmp = xcur;
        xcur = xnew;
		xnew = xtmp;
    }
   
    return xnew;
}

Vector Sparse::conjGradient(Vector const& b) const
{
	assert(this->sizeColumns() == b.size() && this->sizeRows() == b.size());
	
	uint n = sizeColumns();

	Vector x(n);
	Vector r = b; // r = b - A * 0
	Vector p = r;
	
	double rSquare = dot(r,r);
	
	// in the worst case we have to calculate all of the conjugate directions
	for(uint i = 0; i < n; i++)
	{
		double alpha = rSquare / dot(p, (*this) * p);
		Vector v = alpha * p;
		x += v;
		r = r - alpha * ((*this) * p);
		double rSquareNew = dot(r,r);
		if (rSquareNew < CONJGRADIENT_TOLERANCE)
		{
			break;
		}
		p = r + (rSquareNew / rSquare) * p;
		rSquare = rSquareNew;
	}
	
    return x;
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
    int row_cur, row_prev=-1;
    
    for(unsigned int k = 0; k < sizeNNZ; k++)
    {
        is >> row_cur >> colInd[k] >> vals[k];
		if (row_cur > row_prev)
        {
			rowPtr[row_cur] = k;
		}
		row_prev = row_cur;
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