#ifndef VOLUME_H
#define VOLUME_H

#include "utility.h"

using namespace std;
using namespace igl::opengl::glfw;
using namespace t_mesh;
class Volume {
public:
	Volume():id(-1){}
	void setViewer(Viewer* viewer) {
		this->viewer = viewer; // 在头文件中内联，若在cpp中用inline，会导致其他文件找不到定义
	}
	void draw(bool tmesh, bool polygon, bool surface, double resolution = 0.01);
	virtual Point3d eval(double u, double v, double w) = 0;

	virtual int readVolume(string) = 0;
	virtual int saveVolume(string) = 0;
	int saveAsHex(string, double resolution = 0.01);
protected:
	virtual void drawTmesh() = 0;
	virtual void drawControlpolygon() = 0;
	virtual void drawVolume(double resolution = 0.01);

protected:
	Viewer* viewer;
	int id;

};


#endif // !VOLUME_H

