#ifndef TSPLINEVOLUME_H
#define TSPLINEVOLUME_H
#include "utility.h"

using namespace std;
using namespace igl::opengl::glfw;
using namespace t_mesh;
class TsplineVolume {
public:
	void setViewer(Viewer* viewer);
	void draw(bool tmesh, bool polygon, bool surface, double resolution);
	Point3d eval(double u, double v, double w);

	int readVolume(string);
	int saveVolume(string);
private:
	void drawTmesh();
	void drawControlpolygon();
	void drawVolume(double resolution = 0.01);

private:
	map<double, Mesh3d*>       w_map;   
	Eigen::VectorXd            w_knots;     // w向节点向量
	
	Viewer* viewer;
	
};

#endif // !TSPLINEVOLUME_H

