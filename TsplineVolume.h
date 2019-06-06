#ifndef TSPLINEVOLUME_H
#define TSPLINEVOLUME_H
#include "Volume.h"
class TsplineVolume : public Volume{
private:
	TsplineVolume& operator=(const TsplineVolume&){}
public:
	TsplineVolume(){}
	TsplineVolume(const TsplineVolume& other);
	~TsplineVolume();
	Point3d eval(double u, double v, double w) override;
	int readVolume(string) override;
	int saveVolume(string) override;

	void insert(double w, Mesh3d* mesh);

private:
	void drawTmesh() override;
	void drawControlpolygon() override;

public:
	map<double, Mesh3d*>       w_map;   
	Eigen::VectorXd            w_knots;     // w向节点向量

	
};

#endif // !TSPLINEVOLUME_H

