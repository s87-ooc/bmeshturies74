#ifndef __TYPES_H__
#define __TYPES_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: types.h

 Description: définition des types utilisés partout

**************************************************************/

// general
typedef unsigned int uint;

#define EQ_TOL 0.0000000001

// ----------------------------------------------------------------------------

/** get the position in an array based on an elements address */
template <typename T>
uint getIndex(T* start, const T& entry)
{
	return &entry - start;
}

/** delete an array after deleting the objects pointed to */
template <typename T>
void deletePtrArray(T* array, uint size)
{
	if (!array)
	{
		return;
	}
	
	for (uint i = 0; i < size; i++)
	{
		delete array[i];
	}
	
	delete[] array;
}

// ----------------------------------------------------------------------------

/** macro for safe deletion of pointers */
#define SAFE_DELETE(P) if (P) delete P

/** macro for safe deletion of dynamic arrays */
#define SAFE_ARRDELETE(A) if (A) delete[] A

/** macro for dumping an array to the console */
#define DUMP_ARR(A, NUM) cout << #A << ": [ "; for (uint i = 0; i < NUM; i++) { cout << A[i] << " "; } cout << " ]" << endl 

/** macro for dumping a vector to the console */
#define DUMP_VEC(V) cout << #V << ": [ "; for (uint i = 0; i < V.size(); i++) { cout << V(i); if(i == V.size() - 1) cout << " ] "; else cout << " ; "; } cout << V.size() << endl

/** macro for dumping a matrix to the console */
#define DUMP_MAT(M) cout << #M << ": [ "; for (uint i = 0; i < M.sizeRows(); i++) { for (uint j = 0; j < M.sizeColumns(); j++) {	cout << M(i, j) << " ";	} if (i < M.sizeRows() - 1)	{ cout << "; ";	} } cout << "] " << M.sizeRows() << "x" << M.sizeColumns() << endl

#endif // __TYPES_H__
