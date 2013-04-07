/*************************************************************

 Mini Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: visualization.cpp

 Description: implementation des classes pour la visualization

**************************************************************/

#include <iostream>

#include "visualization.h"

using namespace std;

// ----------------------------------------------------------------------------

CPlot::CPlot()
{
	cout << "Plot created" << endl;
}

CPlot::~CPlot()
{
	cout << "Plot destroyed" << endl;
}