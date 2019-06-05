#ifndef PIANASRIMETHOD_H
#define PIANASRIMETHOD_H

#include "PiaMethod.h"
#include <igl/point_mesh_squared_distance.h>
class PiaNasriMethod :public PiaMethod {
public:
	PiaNasriMethod(const vector<NURBSCurve>& _curves, int _maxIterNum = 100, double _eps = 1e-5)
		:PiaMethod(_curves, _maxIterNum, _eps) {

	}
	// 设置辅助点
	void set_helper_points(const MatrixXd& points) {
		helper_points.resize(points.rows());
		for (int i = 0; i < points.rows(); i++) {
			helper_points[i].origin.fromVectorXd(points.row(i));
		}
	}
	// 将辅助点参数化
	void param_helper_points();

	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 按一定规则在误差大的地方插入节点，局部加细
	void calculate() override;  // 计算流程

public:
	void sample_fitPoints() override;
	//virtual void fit();
	void pia() override;

private:
	vector<FitPoint> helper_points;
};
#endif // !PIANASRIMETHOD_H

