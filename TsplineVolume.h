#ifndef TSPLINEVOLUME_H
#define TSPLINEVOLUME_H
#include "utility.h"

using namespace std;
using namespace igl::opengl::glfw;
using namespace t_mesh;
class TsplineVolume {
public:
	void setViewer(Viewer* viewer) {
		this->viewer = viewer; // 在头文件中内联，若在cpp中用inline，会导致其他文件找不到定义
	}
	void draw(bool tmesh, bool polygon, bool surface, double resolution = 0.01);
	Point3d eval(double u, double v, double w);

	int readVolume(string);
	int saveVolume(string);
	int saveAsHex(string, double resolution = 0.01);
private:
	void drawTmesh();
	void drawControlpolygon();
	void drawVolume(double resolution = 0.01);

private:
	map<double, Mesh3d*>       w_map;   
	Eigen::VectorXd            w_knots;     // w向节点向量
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXi Hex;
	Viewer* viewer;
	
};

#endif // !TSPLINEVOLUME_H

