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

// ----------------------------------------------------------------------------

CVector::CVector()
{	
	m_Entries = 0;
	m_Size = 0;
}

CVector::~CVector()
{	
	if (m_Entries > 0)
	{
		delete[] m_Entries;
	}
}

CVector::CVector(uint size, double* data)
{	
	assert (size > 0);
	
	m_Size = size;
	m_Entries = new double[m_Size];
	
	if (data)
	{
		// we're assuming that the data array has a proper size
		for (uint i = 0; i < m_Size; i++)
		{
			m_Entries[i] = data[i];
		}
	}
}

CVector::CVector(const CVector& v)
{
	m_Size = v.size();
	
	m_Entries = new double[m_Size];
	
	for (uint i = 0; i < m_Size; i++)
	{
		m_Entries[i] = v[i];
	}
}

uint CVector::size() const
{
	return m_Size;
}

void CVector::zero()
{
	memset(m_Entries, 0, sizeof(double) * m_Size);
}

const double& CVector::operator[] (uint i) const
{
	assert(i < m_Size);
	
	return m_Entries[i];
}

double& CVector::operator[] (uint i)
{
	assert(i < m_Size);
	
	return m_Entries[i];
}

CVector CVector::operator- ()
{
	CVector v = *this;
	
	for (uint i = 0; i < m_Size; i++)
	{
		v[i] *= -1.0;
	}
	
	return v;
}

CVector& CVector::operator+= (const CVector& v)
{
	assert(v.size() == m_Size);

	for (uint i = 0; i < m_Size; i++)
	{
		m_Entries[i] += v[i];
	}
	
	return *this;
}

CVector& CVector::operator-= (const CVector& v)
{
	assert(v.size() == m_Size);

	for (uint i = 0; i < m_Size; i++)
	{
		m_Entries[i] -= v[i];
	}
	
	return *this;
}


CVector CVector::operator+ (const CVector& v)
{
	assert(v.size() == m_Size);

	CVector vAdd = v;
	
	for (uint i = 0; i < m_Size; i++)
	{
		vAdd[i] += m_Entries[i];
	}
	
	return vAdd;
}

CVector CVector::operator- (const CVector& v)
{
	assert(v.size() == m_Size);

	CVector vSub = *this;
	
	for (uint i = 0; i < m_Size; i++)
	{
		vSub[i] -= v[i];
	}
	
	return vSub;
}

CVector& CVector::operator*= (const double s)
{
	for (uint i = 0; i < m_Size; i++)
	{
		m_Entries[i] *= s;
	}
	
	return *this;
}

CVector operator* (const double s, const CVector& v)
{
	CVector vScaled = v;
	for (uint i = 0; i < vScaled.size(); i++)
	{
		vScaled[i] *= s;
	}
	
	return vScaled;
}

CVector operator* (const CVector& v, const double s)
{
	CVector vScaled = v;
	for (uint i = 0; i < vScaled.size(); i++)
	{
		vScaled[i] *= s;
	}
	
	return vScaled;
}

ostream& operator<<(ostream& o, const CVector& v)
{
	o << "[ ";
    
	for (uint i = 0; i < v.size(); i++)
	{
		o << v[i] << " ";
	}
	
	o << "] " << v.size();
    return o;
} 

/*istream& operator>>(istream& i, CVector& v)
{
    float re,im;
    i >> re >> im;
    c = complexe(re,im);
    return i;
}*/

// ----------------------------------------------------------------------------

const double CMatrixData::s_Zero = 0.0;

CMatrixDataDense::CMatrixDataDense()
{
	m_SizeRows = 0;
	m_SizeColumns = 0;
	m_Data = 0;
}

CMatrixDataDense::~CMatrixDataDense()
{
	if (m_Data)
	{
		delete[] m_Data;
	}
}

CMatrixDataDense::CMatrixDataDense(uint rows, uint columns, double* data)
{
	m_SizeRows = rows;
	m_SizeColumns = columns;
	m_Data = new double[m_SizeRows * m_SizeColumns];
	
	if (data)
	{
		for (uint i = 0; i < m_SizeRows * m_SizeColumns; i++)
		{
			m_Data[i] = data[i];
		}
	}
}

CMatrixDataDense::CMatrixDataDense(const CMatrixDataDense& mdata)
{
	m_SizeRows = mdata.m_SizeRows;
	m_SizeColumns = mdata.m_SizeColumns;
	m_Data = new double[m_SizeRows * m_SizeColumns];
	
	for (uint i = 0; i < m_SizeRows * m_SizeColumns; i++)
	{
		m_Data[i] = mdata.m_Data[i];
	}
}

/*double& CMatrixDataDense::operator() (uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "nonconst ";
	
	return m_Data[m_SizeColumns * r + c];
}

const double& CMatrixDataDense::operator() (uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "const as rValue" << endl;
	
	return m_Data[m_SizeColumns * r + c];
}*/

// @ workaround
double& CMatrixDataDense::_l(uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	return m_Data[m_SizeColumns * r + c];
}

const double& CMatrixDataDense::_r(uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	return m_Data[m_SizeColumns * r + c];
}
// @@@

// ----------------------------------------------------------------------------

CMatrixDataDiag::CMatrixDataDiag()
{
	m_SizeRows = 0;
	m_SizeColumns = 0;
	m_Data = 0;
}

CMatrixDataDiag::~CMatrixDataDiag()
{
	if (m_Data)
	{
		delete[] m_Data;
	}
}

CMatrixDataDiag::CMatrixDataDiag(uint size, double* data)
{
	m_SizeRows = m_SizeColumns = size;
	m_Data = new double[size];
	
	if (data)
	{
		for (uint i = 0; i < m_SizeRows; i++)
		{
			m_Data[i] = data[i];
		}
	}
}

CMatrixDataDiag::CMatrixDataDiag(const CMatrixDataDiag& mdata)
{
	m_SizeRows = m_SizeColumns = mdata.m_SizeRows;
	m_Data = new double[m_SizeRows];
	
	for (uint i = 0; i < m_SizeRows; i++)
	{
		m_Data[i] = mdata.m_Data[i];
	}
}

/*double& CMatrixDataDiag::operator() (uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "nonconst ";

	return m_Data[r];
}

const double& CMatrixDataDiag::operator() (uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "const as rValue" << endl;
	
	if (r == c)
	{
		return m_Data[m_SizeColumns * r + c];
	}
}*/

// @ workaround
double& CMatrixDataDiag::_l(uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);

	if (r != c)
	{
		cout << "lValue only possible on diagonal!" <<endl;
		assert(false);
	};
	
	return m_Data[r];
}

const double& CMatrixDataDiag::_r(uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	if (r != c)
	{
		return CMatrixData::s_Zero;
	}
	else
	{
		return m_Data[r];
	}
}
// @@@

// ----------------------------------------------------------------------------

CMatrixDataSparse::CMatrixDataSparse()
{
	m_SizeRows = 0;
	m_SizeColumns = 0;
	m_Data = 0;
}

CMatrixDataSparse::~CMatrixDataSparse()
{
	if (m_Data)
	{
		delete[] m_Data;
	}
}

CMatrixDataSparse::CMatrixDataSparse(uint rows, uint columns)
{
	m_SizeRows = rows;
	m_SizeColumns = columns;
	
	m_Data = 0;
	// TODO
}

CMatrixDataSparse::CMatrixDataSparse(const CMatrixDataSparse& mdata)
{
	m_SizeRows = mdata.m_SizeRows;
	m_SizeColumns = mdata.m_SizeColumns;
	
	m_Data = 0;
	// TODO
}

/*double& CMatrixDataSparse::operator() (uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "nonconst ";
	
	return m_Data[m_SizeColumns * r + c];
}

const double& CMatrixDataSparse::operator() (uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	//cout << "const as rValue" << endl;
	
	return m_Data[m_SizeColumns * r + c];
}*/

// @ workaround
double& CMatrixDataSparse::_l(uint r, uint c)
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	return m_Data[0];
}

const double& CMatrixDataSparse::_r(uint r, uint c) const
{
	assert (r < m_SizeRows);
	assert (c < m_SizeColumns);
	
	return m_Data[0];
}
// @@@

// ----------------------------------------------------------------------------

CMatrix::CMatrix()
{
	m_Data = 0;
	m_Type = eMDT_DENSE;
}

CMatrix::~CMatrix()
{
	if (m_Data)
	{
		delete m_Data;
	}
}

CMatrix::CMatrix(uint rows, uint columns, eMatrixDataType type, double* data)
{
	assert(rows > 0);
	assert(columns > 0);
	
	m_Type = type;
	
	if (m_Type == eMDT_DIAG && rows != columns)
	{
		cout << "Diagonal matrix has to be a square matrix!" << endl;
		assert(false);
	}
	
	switch(m_Type)
	{
		case eMDT_DENSE:
			m_Data = new CMatrixDataDense(rows, columns, data);
			break;
		case eMDT_DIAG:
			m_Data = new CMatrixDataDiag(rows, data);
			break;
		case eMDT_SPARSE:
			m_Data = new CMatrixDataSparse(rows, columns);
			break;
		default:
			cout << "Invalid matrix type!" << endl;
			assert(false);
			break;
	}
}

CMatrix::CMatrix(const CMatrix& m)
{
	m_Type = m.m_Type;

	if (m_Type == eMDT_DENSE)
	{
		CMatrixDataDense* pData = dynamic_cast<CMatrixDataDense*>(m.m_Data);
		m_Data = new CMatrixDataDense(*pData);
	}
	else if (m_Type == eMDT_DIAG)
	{
		CMatrixDataDiag* pData = dynamic_cast<CMatrixDataDiag*>(m.m_Data);
		m_Data = new CMatrixDataDiag(*pData);
	}
	else if (m_Type == eMDT_SPARSE)
	{
		CMatrixDataSparse* pData = dynamic_cast<CMatrixDataSparse*>(m.m_Data);
		m_Data = new CMatrixDataSparse(*pData);
	}
	else
	{
		cout << "Invalid matrix type!" << endl;
		assert(false);
	}
}

uint CMatrix::sizeRows() const
{
	assert(m_Data);
	
	return m_Data->m_SizeRows;
}

uint CMatrix::sizeColumns() const
{
	assert(m_Data);
	
	return m_Data->m_SizeColumns;
}

/*double& CMatrix::operator() (uint r, uint c)
{
	assert(m_Data);
	
	return (*m_Data)(r, c);
}

const double& CMatrix::operator() (uint r, uint c) const
{
	assert(m_Data);
	
	return (*m_Data)(r, c);
}*/

// @ workaround
double& CMatrix::_l(uint r, uint c)
{
	assert(m_Data);
	
	return m_Data->_l(r, c);
}

const double& CMatrix::_r(uint r, uint c) const
{
	assert(m_Data);
	
	return m_Data->_r(r, c);
}
// @@@

CMatrix CMatrix::operator- ()
{
	assert (m_Data);
	
	CMatrix m(*this);
	
	for (uint i = 0; i < m.sizeRows(); i++)
	{
		for (uint j = 0; j < m.sizeColumns(); j++)
		{
			m._l(i, j) *= -1.0;
		}
	}
	
	return m;
}

CMatrix& CMatrix::operator+= (const CMatrix& m)
{
	assert (m_Data->m_SizeRows == m.sizeRows());
	assert (m_Data->m_SizeColumns == m.sizeColumns());
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			m_Data->_l(i, j) += m.m_Data->_r(i, j);
		}
	}
	
	return *this;
}

CMatrix& CMatrix::operator-= (const CMatrix& m)
{
	assert (m_Data->m_SizeRows == m.sizeRows());
	assert (m_Data->m_SizeColumns == m.sizeColumns());
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			m_Data->_l(i, j) -= m.m_Data->_r(i, j);
		}
	}
	
	return *this;
}

CMatrix CMatrix::operator+ (const CMatrix& m)
{
	assert (m_Data->m_SizeRows == m.sizeRows());
	assert (m_Data->m_SizeColumns == m.sizeColumns());
	
	CMatrix mAdd(*this);
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			mAdd.m_Data->_l(i, j) += m.m_Data->_r(i, j);
		}
	}
	
	return mAdd;
}

CMatrix CMatrix::operator- (const CMatrix& m)
{
	assert (m_Data->m_SizeRows == m.sizeRows());
	assert (m_Data->m_SizeColumns == m.sizeColumns());
	
	CMatrix mSub(*this);
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			mSub.m_Data->_l(i, j) -= m.m_Data->_r(i, j);
		}
	}
	
	return mSub;
}

CMatrix& CMatrix::operator*= (const double s)
{
	assert (m_Data);
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			m_Data->_l(i, j) *= s;
		}
	}
	
	return *this;
}

/*CMatrix operator* (const CVector& v, const CMatrix& m)
{
	assert(v.size() == m.sizeRows());

	CMatrix mMult(m);
	
	for (uint i = 0; i < m_Data->m_SizeRows; i++)
	{
		for (uint j = 0; j < m_Data->m_SizeColumns; j++)
		{
			m.m_Data->_l(i, j) = 0.0;
			for (uint k = 0; k < v.size(); k++)
			{
				m.m_Data->_l(i, j) += v[k] * m.m_Data->_l(k, j);
			}
		}
	}
	
	return mMult;
}*/

/*
CVector operator* (const CMatrix& m, const CVector& v);
CMatrix operator* (const double s, const CMatrix& m);
CMatrix operator* (const CMatrix& m, const double s);
*/

ostream& operator<<(ostream& o, const CMatrix& m)
{
	o << "[ ";
    
	for (uint i = 0; i < m.sizeRows(); i++)
	{
		for (uint j = 0; j < m.sizeColumns(); j++)
		{
			//o << m(i, j) << " ";
			o << m._r(i, j) << " ";
		}
		
		if (i < m.sizeRows() - 1)
		{
			cout << "; ";
		}
	}
	
	o << "] " << m.sizeRows() << "x" << m.sizeColumns() << " ";
    
	return o;
}