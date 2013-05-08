/*************************************************************

 Mini Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
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

Vector& Vector::constructFunc(const Mesh& mesh, double (*f)(const Vertex&))
{
	assert(size() == mesh.countVertices());
	
	for (uint i = 0; i < size(); i++)
	{
		mVals[i] = f(mesh.V[i]);
	}
	
	return *this;
}

Vector& Vector::constructFuncIntSurf(const Mesh& mesh, double (*f)(const Vertex&))
{
	assert(size() == mesh.countVertices());
	
	for (uint e = 0; e < mesh.countEdges(); e++)
	{		
		// interpolation by affine function, only the function contributes
		double factor = (mesh.E[e]).length * 0.5;
		
		mVals[mesh.E[e](0).id] += factor * f(mesh.E[e](0));
		mVals[mesh.E[e](1).id] += factor * f(mesh.E[e](1));
		
		//cout << "e: " << e << " " << mesh.E[e](0).id << " " << mesh.E[e](1).id << " " << factor * f(mesh.E[e](0)) << " " << mesh.E[e].length << endl;
	}
	
	return *this;
}

istream& operator>> (istream &in, Vector& v)
{
	v.mVals.clear();
	
	uint size;
	double val;
	
	in >> size;
	
	for (uint i = 0; i < size; i++)
	{
		in >> val;
		v.mVals.push_back(val);
	}
	
	return in;
}

ostream& operator<<(ostream &out, Vector& v)
{
	out << v.size() << endl;
	for (uint i = 0; i < v.size() - 1; i++)
	{
		out << v(i) << " ";
	}
	out << v(v.size() - 1) << endl;
	
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

// double& SparseMap::operator() (uint row, uint col)
// {
// 	return mVals[row*mSizeColumns + col];
// }

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
	assert(rhs.sizeColumns() == this->sizeColumns());
	assert(rhs.sizeRows() == this->sizeRows());

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

SparseMap& SparseMap::constructA(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	uint Nt = mesh.countTriangles();

	TTriangles const& T = mesh.T;

	mSizeRows = Nv;
	mSizeColumns = Nv;

	/*double a,b,c,d, A, B, C, D, l11, l22, l33, l12, l13, l23;
	// iterate over the triangles of the mesh
	for( uint t = 0; t < Nt; t++)
	{
		a = T[t](1).x - T[t](0).x;
		b = T[t](2).x - T[t](0).x;
		c = T[t](1).y - T[t](0).y;
		d = T[t](2).y - T[t](0).y;
		A = pow(a,2) + pow(c,2);
		B = pow(b,2) + pow(d,2);
		C = a*b + c*d;
		D = fabs(a*d - b*c);

		l11 = (A + B - 2*C) / D;
		l22 = B / D;
		l33 = A / D;
		l23 = - (C / D);
		l13 = (C - A) / D;
		l12 = (C - B) / D;

		this->addAt(T[t](0).id, T[t](0).id, l11);
		this->addAt(T[t](1).id, T[t](1).id, l22);
		this->addAt(T[t](2).id, T[t](2).id, l33);
		this->addAt(T[t](0).id, T[t](1).id, l12);
		this->addAt(T[t](1).id, T[t](0).id, l12);
		this->addAt(T[t](0).id, T[t](2).id, l13);
		this->addAt(T[t](2).id, T[t](0).id, l13);
		this->addAt(T[t](1).id, T[t](2).id, l23);
		this->addAt(T[t](2).id, T[t](1).id, l23);
	}
	*/
	double x21, x31, y21, y31, a, b, d, D, l11, l22, l33, l12, l13, l23;
	// iterate over the triangles of the mesh
	for( uint t = 0; t < Nt; t++)
	{
		//cout << t << endl;
	
		x21 = T[t](1).x - T[t](0).x;
		x31 = T[t](2).x - T[t](0).x;
		y21 = T[t](1).y - T[t](0).y;
		y31 = T[t](2).y - T[t](0).y;
		
		D = fabs(x21*y31 - x31*y21);
		
		a = (pow(x31, 2) + pow(y31, 2)) / D;
		d = (pow(x21, 2) + pow(y21, 2)) / D;
		b = - (x31 * x21 + y21 * y31) / D;
		
		l11 = a + d + 2*b;
		l22 = a;
		l33 = d;
		l12 = -a - b;
		l13 = -b - d;
		l23 = b;
		
		addAt(T[t](0).id, T[t](0).id, l11);
		addAt(T[t](1).id, T[t](1).id, l22);
		addAt(T[t](2).id, T[t](2).id, l33);
		addAt(T[t](0).id, T[t](1).id, l12);
		addAt(T[t](1).id, T[t](0).id, l12);
		addAt(T[t](0).id, T[t](2).id, l13);
		addAt(T[t](2).id, T[t](0).id, l13);
		addAt(T[t](1).id, T[t](2).id, l23);
		addAt(T[t](2).id, T[t](1).id, l23);
	}

	return (*this);
}

SparseMap& SparseMap::constructM(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	uint Nt = mesh.countTriangles();

	TTriangles const& T = mesh.T;

	mSizeRows = Nv;
	mSizeColumns = Nv;

	double mel = 1/12.0; // since det = 2 * area we have 2/24 = 1/12
	double value;
	double coef;

	// iterate over the triangles of the mesh
	for( uint t = 0; t < Nt; t++)
	{

		// iterate over vertices of the current triangle
		for(uint i=0; i < 3; i++)
		{
			for(uint j=0; j < 3; j++)
			{
				coef = ( i == j ) ? 2. : 1.;
				value = coef * mel * T[t].area;
				addAt(T[t](i).id, T[t](j).id, value);
			}
		}
	}

	return (*this);
}


SparseMap& SparseMap::constructB(const Mesh& mesh)
{
	uint Nv = mesh.countVertices();
	uint Ne = mesh.countEdges();

	TEdges const& E = mesh.E;

	mSizeRows = Nv;
	mSizeColumns = Nv;

	double value;

	for( uint e = 0; e < Ne; e++)
	{
		value = E[e].length / 6.0;
		/*this->addAt(E[e](0).id, E[e](0).id, value);
		this->addAt(E[e](1).id, E[e](1).id, value);
		this->addAt(E[e](0).id, E[e](1).id, value);
		this->addAt(E[e](1).id, E[e](0).id, value);*/
		// stjepan: integral for i = j gives E[e].length / 3.0
		addAt(E[e](0).id, E[e](0).id, value * 2.);
		addAt(E[e](1).id, E[e](1).id, value * 2.);
		addAt(E[e](0).id, E[e](1).id, value);
		addAt(E[e](1).id, E[e](0).id, value);
	}

	return (*this);
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

SparseLIL SparseLIL::prodTranspose() const
{
	assert (sizeRows() == sizeColumns());
	
	SparseLIL m(sizeRows(), sizeColumns());
	
	//val += (i, k) * (j, k)
	
	for (uint i = 0; i < sizeRows(); i++)
	{
		if (i % 10 == 0)
			cout << "i " << i << endl;
		
		for (uint j = 0; j < sizeColumns(); j++)
		{
			//if (j % 100 == 0)
			//	cout << "j " << j << endl;
			// TODO: this may be VERY inefficient!
			double val = 0.;
			for (uint k = 0; k < sizeColumns(); k++)
			{
				for (uint l = 0; mColInd[j][l] <= k; l++)
				{
					//val += (*this)(i, k) * (*this)(j, k);
					if (mColInd[j][l] == k)
					{
						for (uint m = 0; mColInd[i][m] <= k; m++)
						{
							if (mColInd[i][m] == k)
							{
								val += (*this)(i, k) * mRowVals[j][l];
							}
						}
					}
				}
			}
			if (fabs(val) >= EQ_TOL)
			{
				m(i, j) = val;
			}
		}
	}
	
	return m;
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
	int index, row_cur, row_prev = -1;
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

	// last row points to number of nonzero values (last nnz index + 1, as starts with 0)
	mRowPtr[M.sizeRows()] = M.sizeNNZ();
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

Vector Sparse::LU(Vector const& b, Sparse* lu) const
{
	const uint n = sizeRows();

	// store the permutations of partial pivoting
	uint* p = new uint[n];
	for (uint i = 0; i < n; i++)
	{
		p[i] = i;
	}

	SparseLIL A(*this);
	
	// ----------

	// decompose A = LU with partial pivoting
	
	for (uint k = 0; k < n - 1; k++)
	{
		//cout << "-" << k << " " << A(p[k], p[k]) << endl;
	
		/*cout << "--- STEP ---" << endl;
		for (uint i = 0; i < n; i++)
		{
			for (uint j = 0; j < n; j++)
			{
				cout << A(p[i], j) << " ";
			}
			cout << endl;
		}*/
		
		if (fabs(A(p[k], k)) < EQ_TOL)
		{
			//cout << "ZERO " << A(p[k], k) << endl;
			// TODO: perform permutation
			for (uint k0 = k+1; k0 < n; k0++)
			{
				//cout << k0 << " " << A(p[k0], k) << endl;
				if (fabs(A(p[k0], k)) >= EQ_TOL)
				{
					//cout << p[k] << " <-> " << p[k0] << endl;
					uint swap = p[k];
					p[k] = k0;
					p[k0] = swap;
					break;
				}
			}
		}

		/*cout << "--- SWAPPED ---" << endl;
		for (uint i = 0; i < n; i++)
		{
			for (uint j = 0; j < n; j++)
			{
				cout << A(p[i], j) << " ";
			}
			cout << endl;
		}*/
		
		for (uint i = k+1; i < n; i++)
		{
			A(p[i],k) = A(p[i],k) / A(p[k],k);
			
			for (uint j = k+1; j < n; j++)
			{
				A(p[i],j) = A(p[i],j) - A(p[i],k)*A(p[k],j);
			}
		}
	}
	
	// save the decomposition in m if the pointer was valid
	if (lu)
	{
		*lu = Sparse(A);
	}
	
	// ----------
	
	// Reconstruct solution out of LU decomposition
	
	Vector x(n);
	Vector y(n);

	y(p[0]) = b(p[0]);	// do this outside of the loop as we want to go over uint, -1 not uint
	for (uint k = 1; k < n; k++)
	{
		y(p[k]) = b(p[k]);
		for (uint j = 0; j <= k-1; j++)
		{
			y(p[k]) -= A(p[k], j) * y(p[j]);
		}
	}
	
	//cout << "y: " << y << endl;
	
	for (uint k = n-1; k >= 1; k--)
	{
		//cout << k << " " << k+1 << " " << n-1 << endl;
		x(k) = y(p[k]);
		for (uint j = n-1; j >= k+1; j--)
		{
			//cout << "GO" << endl;
			x(k) -= A(p[k], j) * x(j);
		}
		x(k) /= A(p[k], k);
	}
	// do this outside of the loop since we're using uint
	//cout << 0 << " " << 0+1 << " " << n-1 << endl;
	x(0) = y(p[0]);
	for (uint j = n-1; j >= 1; j--)
	{
		//cout << "GO" << endl;
		x(0) -= A(p[0], j) * x(j);
	}
	x(0) /= A(p[0], 0);
	
	//cout << "x: " << x << endl;
	
	// ----------
	
	// cleanup memory
	delete[] p;	
	
	return x;
}

// ----------------------------------------------------------------------------

istream& operator>> (istream &is, Sparse &m)
{
	unsigned int sizeRows, sizeColumns, sizeNNZ;
	is >> sizeRows >> sizeColumns >> sizeNNZ;

    TVals vals(sizeNNZ);    
	TInd colInd(sizeNNZ);
	TInd rowPtr(sizeRows + 1);
    
    /*unsigned int r;
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
    
    rowPtr[sizeRows] = rowPtr[sizeRows - 1] + 1;*/
	
	// stjepan: use compact format for sparse matrices
	
	for (uint i = 0; i < vals.size(); i++)
	{
		is >> vals[i];
	}
	
	for (uint i = 0; i < colInd.size(); i++)
	{
		is >> colInd[i];
	}
	
	for (uint i = 0; i < rowPtr.size() - 1; i++)
	{
		is >> rowPtr[i];
	}
	rowPtr[rowPtr.size() - 1] = vals.size();
    
    m = Sparse(vals, colInd, rowPtr, sizeColumns);
    
    return is;
}

ostream& operator<< (ostream &os, Sparse &m)
{
	os << m.sizeRows() << " "
		<< m.sizeColumns() << " "
		<< m.sizeNNZ() << endl;
	
	/*for (unsigned int i = 0; i < m.sizeRows(); i++)
	{
		for (unsigned int k = m.mRowPtr[i]; k < m.mRowPtr[i+1]; k++)
		{
			os << i << " " << m.mColInd[k] << " " << m.mVals[k] << endl;
		}
	}*/
	
	// stjepan: store sparse matrices in a compact manner, use DUMP_MAT for printing the values in a human-readable manner
	
	for (uint i = 0; i < m.mVals.size() - 1; i++)
	{
		os << m.mVals[i] << " ";
	}
	os << m.mVals[m.mVals.size() - 1] << endl;
	
	for (uint i = 0; i < m.mColInd.size() - 1; i++)
	{
		os << m.mColInd[i] << " ";
	}
	os << m.mColInd[m.mColInd.size() - 1] << endl;
	
	for (uint i = 0; i < m.mRowPtr.size() - 2; i++)
	{
		os << m.mRowPtr[i] << " ";
	}
	os << m.mRowPtr[m.mRowPtr.size() - 2] << endl;

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