#ifndef VOLUMEPIAMETHOD_H
#define VOLUMEPIAMETHOD_H

#include "VolumeSkinning.h"
#include "FitPoint.hpp"
#include "BsplineVolume.h"
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
	void sample_fitPoints_bvolume();
	// 生成形状引导线，分段优化生成Jacobian值全为正的B样条体，采样用于拟合的点
	void sample_fitPoints_multiVolume();
	void sample_fitPoints();
	void fit();
	void pia();
	void cal_basis_cache();

private:
	const int maxIterNum;
	const double eps;
	double error;
	vector<FitPoint3D> fitPoints;
	vector<FitPoint3D> surface_points;
	vector<FitPoint3D> inter_points;
	vector<FitPoint3D> helper_points;
	vector<vector<double>> basis_cache; // (B_i(t_j)
	vector<double> basis_cache_sum; // sum_j (B_i(t_j)
};
#endif // !VOLUMEPIAMETHOD_H

