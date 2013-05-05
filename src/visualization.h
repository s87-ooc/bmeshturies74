#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: visualization.h

 Description: déclaration des classes pour la visualization

**************************************************************/

// ----------------------------------------------------------------------------

/** which application is used for plotting */
enum EPlotType
{
	ePT_GNUPLOT = 0,	/** Gnuplot: (v.x, v.y, value in v) */
	ePT_MEDIT			/** Medit: (v.x, v.y, value in v) */
};

enum EPlotData
{
	ePD_MESH = 0,	/** no extra data, we're only interested in plotting the mesh */
	ePD_FUNCTION,	/** function pointer is specified that gets evaluated in vertices */
	ePD_VECTOR		/** vector with values for each vertex is specified */
};

// ----------------------------------------------------------------------------

/** Plot a 2-dimensional graph or a set of points in the plane */
class Plot
{
};

// ----------------------------------------------------------------------------

class Mesh;
class Vertex;
class Vector;

/** Plot a 3-dimensional graph over a mesh or simply the mesh geometry */
class PlotMesh
{
private:
	/** name base for data and script files */
	std::string mName;
	
	/** initial template for the gnuplot script file */
	std::string mTemplate;
	
	/** additional plotting arguments */
	std::string mArgs;
	
	/** mesh the values should be plotted on */
	Mesh* mMesh;
	
	/** type of the plot */
	EPlotType mType;
	
	/** type of the data to plot over a mesh */
	EPlotData mDataType;
	
	/** pointer to values of vertices stored in a vector */
	Vector* mVectorPtr;
	
	/** pointer to a function, will be evaluated in mesh vertices */
	double (*mFuncPtr)(const Vertex&);

public:
	PlotMesh();
	~PlotMesh();
	
	// TODO: simpler plots that work with x, y for errors, etc
	
	/** set up plot based on a mesh only */
	PlotMesh(const char* name, Mesh* msh, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a value vector */
	PlotMesh(const char* name, Mesh* msh, Vector* vals, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a value vector */
	PlotMesh(const char* name, Mesh* msh, double (*func)(const Vertex&), EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** generate the plot out of the mesh and data */
	void generate(bool run = false);
};

#endif // __VISUALIZATION_H__