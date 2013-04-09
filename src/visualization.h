#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@gmail.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: visualization.h

 Description: déclaration des classes pour la visualization

**************************************************************/

// ----------------------------------------------------------------------------

/** which application and plot type is used for plotting */
enum EPlotType
{
	ePT_GNUPLOT_SURFACE = 0,
	ePT_GNUPLOT_ERROR,
	ePT_TECPLOT_SURFACE,
	ePT_TECPLOT_ERROR,
};

// ----------------------------------------------------------------------------

class Mesh;

class Plot
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
	
	/** reference to data that's plotted over the mesh (TODO) */
	// TO BE DEFINED HOW WE STORE THIS / WHAT IS NEEDED

public:
	Plot();
	~Plot();
	
	// TODO: simpler plots that work with x, y for errors, etc
	
	/** plot constructor with meshes */
	Plot(const char* name, const char* temp, Mesh* msh, EPlotType type = ePT_GNUPLOT_SURFACE, const char* args = 0);
	
	/** generates a data and script file depending on type, optionally a viewer is called instantly */
	void generate(bool run = false);
};

#endif // __VISUALIZATION_H__