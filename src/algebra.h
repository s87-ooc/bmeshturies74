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

class Sparse;

class Vector
{
private:
	TVals vals_;
    bool col_vect;
	
public:
	Vector();
    Vector(unsigned int);
	Vector(TVals &vals);
    
	unsigned int size() const;
	double norm2() const;
	
	double operator()(unsigned int pos) const;
	double& operator()(unsigned int pos);
	
	Vector operator*(Sparse const& m) const;
	Vector operator*(double s) const;
	Vector operator+(Vector const& v) const;
	Vector operator-(Vector const& v) const;
	
	Vector& operator+=(Vector& v);
	Vector& operator-=(Vector& v);
};

ostream& operator<<(ostream &out, Vector& v);


class Sparse {
private:
    TVals vals;
    TInd col_ind, row_ptr;
    
    unsigned int nrows_, ncols_, nnz_;
    
    friend istream& operator >>(istream&, Sparse&);
    friend ostream& operator <<(ostream&, Sparse&);
    
public:
    Sparse();
    Sparse(unsigned int nrows, unsigned int ncols);
    
    unsigned int ncols() const;
    unsigned int nrows() const;
    unsigned int nnz() const;
    
    
    double operator()(unsigned int, unsigned int) const;
    double& operator()(unsigned int, unsigned int);
    
    Vector operator*(Vector const&) const;
    
    Vector jacobi(Vector const&) const;
    Vector gradient_conj(Vector const&) const;
    
    
};


#endif // __ALGEBRA_H__