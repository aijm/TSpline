#ifndef MESHRENDER_H
#define MESHRENDER_H


#include "Window.h"
#include "utility.h"
class MeshRender : public Window {
public:
	MeshRender(t_mesh::Mesh3d* _mesh,
		bool _showmesh = false, bool _showpolygon = true,
		bool _showsurface = true, double _resolution = 0.01)
		:mesh(_mesh), showmesh(_showmesh), showpolygon(_showpolygon),
		showsurface(_showsurface), resolution(_resolution)
	{
		mesh->setViewer(&viewer);
	}

protected:
	void draw() override;

private:
	t_mesh::Mesh3d* mesh;
	bool showmesh;
	bool showpolygon;
	bool showsurface;
	double resolution;

};

#endif // !MESHRENDER_H

