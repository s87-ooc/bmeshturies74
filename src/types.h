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

/** */
template <typename T>
uint getIndex(T* start, const T& entry)
{
	return &entry - start;
}

#endif // __TYPES_H__
