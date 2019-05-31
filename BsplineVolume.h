#pragma once
#ifndef BSPLINEVOLUME_H
#define BSPLINEVOLUME_H
#include "Volume.h"
class BsplineVolume :public Volume{
public:

	Point3d eval(double u, double v, double w) override;

	int readVolume(string) override;
	int saveVolume(string) override;

private:
	void drawTmesh() override;
	void drawControlpolygon() override;

private:
	vector<vector<vector<Point3d>>> control_grid;
	vector<Eigen::VectorXd> knot_vector;
	
};
#endif // !BSPLINEVOLUME_H

