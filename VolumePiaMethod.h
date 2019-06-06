#ifndef VOLUMEPIAMETHOD_H
#define VOLUMEPIAMETHOD_H

#include "VolumeSkinning.h"
#include "FitPoint.hpp"
class VolumePiaMethod :public VolumeSkinning {
public:
	VolumePiaMethod(const vector<Mesh3d>& _surfaces,int _maxIterNum=100,double _eps = 1e-5)
		:VolumeSkinning(_surfaces),maxIterNum(_maxIterNum),eps(_eps) {

	}

	void calculate() override;  // 计算流程


public:
	// 设置辅助点
	void set_helper_points(const MatrixXd& points) {
		helper_points.resize(points.rows());
		for (int i = 0; i < points.rows(); i++) {
			helper_points[i].origin.fromVectorXd(points.row(i));
		}
		(*viewer).data().add_points(points, red);
	}
	// 将辅助点参数化
	void param_helper_points(Point3d& low, Point3d& high);
	void sample_fitPoints_2();
	void sample_fitPoints();
	void fit();
	void pia();

private:
	const int maxIterNum;
	const double eps;
	double error;
	vector<FitPoint3D> fitPoints;
	vector<FitPoint3D> surface_points;
	vector<FitPoint3D> inter_points;
	vector<FitPoint3D> helper_points;
};
#endif // !VOLUMEPIAMETHOD_H

