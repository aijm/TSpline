#ifndef PIAMINJAEMETHOD_H
#define PIAMINJAEMETHOD_H

#include "PiaMethod.h"
#include "NURBSSurface.h"
#include <igl/point_mesh_squared_distance.h>
class PiaMinJaeMethod : public PiaMethod {
public:
	PiaMinJaeMethod(const vector<NURBSCurve>& _curves, int _maxIterNum = 100, double _eps = 1e-5)
		:PiaMethod(_curves, _maxIterNum, _eps) {

	}

	// 设置辅助点
	void set_helper_points(const MatrixXd& points) {
		helper_points.resize(points.rows());
		for (int i = 0; i < points.rows(); i++) {
			helper_points[i].origin.fromVectorXd(points.row(i));
		}
		(*viewer).data().add_points(points, red);
	}
	// 将辅助点参数化
	std::tuple<double, double, double, double> param_helper_points();

	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 按一定规则在误差大的地方插入节点，局部加细
	void calculate() override;  // 计算流程

public:
	void sample_fitPoints_1();
	void sample_fitPoints_2();
	void sample_fitPoints() override;
	void pia() override;

private:
	vector<FitPoint> curve_points;
	vector<FitPoint> inter_points;
	vector<FitPoint> helper_points;
};

#endif // !PIAMINJAEMETHOD_H

