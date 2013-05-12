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
	ePT_GNUPLOT = 0,	/** Gnuplot: plot values over the mesh (3D) or a vector (2D) */
	ePT_GNUPLOT_SURF,	/** Gnuplot: plot linear combinations of base functions over a grid on the mesh */
	ePT_MEDIT,				/** Medit: generate a file defining the mesh and the solution (optional) */
	ePT_MATLAB,			/** Matlab: plot values over the mesh (3D) or a vector (2D) */
	ePT_MATLAB_SURF	/** Matlab: plot linear combinations of base functions over a grid on the mesh */
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
	Plot(const char* name, Vector& x, Vector& y, const char* title = 0, const char* templ = 0, const char* args = 0);
	
	/** set up plot based on a function */
	Plot(const char* name, Vector& x, double (&func)(double), const char* title = 0, const char* templ = 0, const char* args = 0);

	/** generate the plot out of the mesh and data */
	void generate(EPlotType type = ePT_GNUPLOT, bool run = false, bool savePNG = false);
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
	PlotMesh(const char* name, Mesh& msh, const char* title = 0);
	
	/** set up plot based on a value vector */
	PlotMesh(const char* name, Mesh& msh, Vector& vals, const char* title = 0);
	
	/** set up plot based on a function */
	PlotMesh(const char* name, Mesh& msh, double (&func)(const Vertex&), const char* title = 0);
	
	/** generate the plot out of the mesh and data in one of the various formats, specify the grid density for surface plots */
	void generate(EPlotType type = ePT_GNUPLOT, bool run = false, bool savePNG = false, const char* templ = 0, const char* args = 0, uint grid = 10);
};

// ----------------------------------------------------------------------------

/** video generation from plotted frames */
enum EPlotVideo
{
	ePV_RAW = 0		/** uncompressed video */
};

namespace PlotVideo
{
	void renderVideo(const char* file, const char* framePrefix, double framerate = 24., const char* title = 0,
		uint width = 640, uint height = 480, EPlotVideo codec = ePV_RAW);
};

#endif // __VISUALIZATION_H__