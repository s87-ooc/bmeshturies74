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
	ePT_GNUPLOT = 0,	/** Gnuplot: (v.x, v.y) for 2d, (v.x, v.y, value in v) for 3d */
	ePT_MEDIT,				/** Medit: */
	ePT_OCTAVE,			/** Octave: */
};

enum EPlotData
{
	ePD_MESH = 0,		/** no extra data, we're only interested in plotting the mesh */
	ePD_FUNCTION,	/** function pointer is specified that gets evaluated in vertices */
	ePD_VECTOR		/** vector with values for each vertex is specified */
};

// ----------------------------------------------------------------------------

/** Plot a 2-dimensional graph or a set of points in the plane */
class Plot
{
private:
	/** name base for data and script files */
	std::string mName;
	
	/** title for the plot */
	std::string mTitle;
	
	/** initial template for the gnuplot script file */
	std::string mTemplate;
	
	/** additional plotting arguments */
	std::string mArgs;
	
	/** type of the plot */
	EPlotType mType;

	// ---
	
	/** type of the data to plot over the x values */
	EPlotData mDataType;
	
	/** x values */
	Vector* mXPtr;
	
	/** y values */
	Vector* mYPtr;
	
	/** ponter to a function to plot */
	double (*mFuncPtr)(double);

public:
	Plot();
	~Plot();
	
	/** set up plot based on a value vector */
	Plot(const char* name, Vector& x, Vector& y, const char* title = 0, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a function */
	Plot(const char* name, Vector& x, double (&func)(double), const char* title = 0, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);

	/** generate the plot out of the mesh and data */
	void generate(bool run = false);
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
	
	/** title for the plot */
	std::string mTitle;
	
	/** initial template for the gnuplot script file */
	std::string mTemplate;
	
	/** additional plotting arguments */
	std::string mArgs;
	
	/** mesh the values should be plotted on */
	Mesh* mMeshPtr;
	
	/** type of the plot */
	EPlotType mType;
	
	// ---
	
	/** type of the data to plot over a mesh */
	EPlotData mDataType;
	
	/** pointer to values of vertices stored in a vector */
	Vector* mVectorPtr;
	
	/** pointer to a function, will be evaluated in mesh vertices */
	double (*mFuncPtr)(const Vertex&);

public:
	PlotMesh();
	~PlotMesh();
	
	/** set up plot based on a mesh only */
	PlotMesh(const char* name, Mesh& msh, const char* title = 0, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a value vector */
	PlotMesh(const char* name, Mesh& msh, Vector& vals, const char* title = 0, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a function */
	PlotMesh(const char* name, Mesh& msh, double (&func)(const Vertex&), const char* title = 0, EPlotType type = ePT_GNUPLOT, const char* templ = 0, const char* args = 0);
	
	/** generate the plot out of the mesh and data */
	void generate(bool run = false);
};

#endif // __VISUALIZATION_H__