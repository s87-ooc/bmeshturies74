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


Vector::Vector(unsigned int size):
vals_(size), col_vect(true)
{
}

Vector::Vector(TVals &vals):
vals_(vals), col_vect(true)
{
}

unsigned int Vector::size() const
{
	return (int)vals_.size();
}


double Vector::norm2() const
{
	double norm = 0;
	
	for (int i = 0; i < vals_.size(); i++)
		norm += pow(vals_[0], 2);
	
	return sqrt(norm);
}

double Vector::operator()(unsigned int pos) const
{
	return vals_[pos];
}

double& Vector::operator()(unsigned int pos)
{
	return vals_[pos];
}

Vector Vector::operator*(Sparse const& M) const
{
	int n = (*this).size();
	Vector res(n);
	
    // implement multiplication
    
	
	return res;
}

Vector Vector::operator*(double s) const
{
    int n = (int)vals_.size();
	Vector res(n);
	
	for (int i = 0; i < n; i++)
		res(i) = s * vals_[i];
	
	return res;
}

Vector Vector::operator+(Vector const& v) const
{
    int n = v.size();
	Vector res(n);
	
	for (int i = 0; i < n; i++)
		res(i) = vals_[i] + v(i);
	
	return res;
}


Vector Vector::operator-(Vector const& v) const
{
    int n = v.size();
	Vector res(n);
	
	for (int i = 0; i < n; i++)
		res(i) = vals_[i] - v(i);
	
	return res;
}

Vector& Vector::operator+=(Vector &v)
{
	for (int i = 0; i < v.size(); i++)
		vals_[i] += v(i);
    
	return *this;
}

Vector& Vector::operator-=(Vector &v)
{
	for (int i = 0; i < v.size(); i++)
		vals_[i] -= v(i);
	
	return *this;
}

// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------


Sparse::Sparse()
{
}

Sparse::Sparse(unsigned int nrows, unsigned int ncols):
nrows_(nrows), ncols_(ncols)
{
}

unsigned int Sparse::ncols() const
{
    return ncols_;
}

unsigned int Sparse::nrows() const
{
    return nrows_;
}

unsigned int Sparse::nnz() const
{
    return nnz_;
}

istream& operator >>(istream &is, Sparse &M)
{
    is >> M.nrows_ >> M.ncols_ >> M.nnz_;
    
    M.col_ind = TInd(M.nnz_);
    M.row_ptr = TInd(M.nrows_+1, -1);
    M.vals = TVals(M.nnz_);
    
    int r;
    
    for(int k = 0; k < M.nnz_; k++)
    {
        is >> r >> M.col_ind[k] >> M.vals[k];
        if (M.row_ptr[r] < 0) { M.row_ptr[r] = k; }
    }
    
    M.row_ptr[M.nnz_] = M.row_ptr[M.nnz_ - 1] + 1;
    
    return is;
}

ostream& operator <<(ostream &os, Sparse &M)
{
    return os;
}


double Sparse::operator()(unsigned int row, unsigned int col) const
{
    cout << "operator() const" << endl;
    if (row > nrows_-1 || col > ncols_-1)
    {
        cout << "Erreur de dimension de la matrice" << endl;
        return 0.0;
    }
    
    for (int i = row_ptr[row]; i < row_ptr[row+1]; i++)
    {
        if ( col_ind[i] == col)
            return vals[i];
    }
    
    return 0.0;
}

double& Sparse::operator()(unsigned int row, unsigned int col)
{
    cout << "operator() non-const" << endl;
    
    if (row > nrows_-1 || col > ncols_-1)
    {
        cout << "Erreur de dimension de la matrice" << endl;
        return vals[0];
    }
    
    
    return vals[0];
}



Vector Sparse::operator*(Vector const& x) const
{
    // dimension check
    // if (ncols != v.size()) { // throw error }
    
    TVals y(ncols_);
    for(int i=0; i<nrows_; i++)
    {
        y[i] = 0.0;
        for(int k = row_ptr[i]; k < row_ptr[i+1]; k++)
        {
            y[i] += vals[k]*x(col_ind[k]);
        }
    }
    return y;
}

Vector Sparse::jacobi(Vector const& b) const
{
    Vector x1(ncols_);
    Vector x2(ncols_);
    Vector& xcur = x1;
    Vector& xnew = x2;
    
    int k = 0;
    double diag;
    bool convergence = false;
    
    while (!convergence)
    {
        for(int i=0; i<nrows_; i++)
        {
            xnew(i) = b(i);
            diag = 0;
            // iterate over the ith row of the matrix
            for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
            {
                if( col_ind[k] != i)
                    xnew(i) -= vals[k]*xcur(col_ind[k]);
                else
                    diag = vals[k];
            }
            // do zero check
            xnew(i) /= diag;
        }
        // TODO: convergence check
        // implement equality operator for Vector
        xcur = xnew;
    }
    
    return xnew;
    
}



Vector Sparse::gradient_conj(Vector const& b) const
{
    Vector p = b - (*this)*b;
    return p;
}
