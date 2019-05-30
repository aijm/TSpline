#ifndef MINJAEMETHOD_H
#define MINJAEMETHOD_H

#include "Skinning.h"
class MinJaeMethod : public Skinning {
public:
	MinJaeMethod(
		const vector<NURBSCurve>& _curves, 
		int _sampleNum = 10, int _maxIterNum = 20, double _eps = 1e-5)
		:Skinning(_curves),sampleNum(_sampleNum),
		maxIterNum(_maxIterNum),eps(_eps) {

	}

	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 在T-preimage中插入中间节点
	void calculate() override;  // 计算流程

private:
	void inter_init();          // 初始化中间点的坐标
	double inter_update();        // 更新中间点的坐标

private:
	const int sampleNum;
	const int maxIterNum;
	const double eps;
	map<double, map<double, Point3d>> initial_cpts;
};


#endif // !MINJAEMETHOD_H