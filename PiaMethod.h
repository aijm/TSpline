#ifndef PIAMETHOD_H
#define PIAMETHOD_H

#include "skinning.h"
#include "FitPoint.hpp"

class PiaMethod : public Skinning {
public:
	PiaMethod(const vector<NURBSCurve>& _curves, int _maxIterNum = 100, double _eps = 1e-5)
		:Skinning(_curves), maxIterNum(_maxIterNum), eps(_eps){

	}
	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 按一定规则在误差大的地方插入节点，局部加细
	void calculate() override;  // 计算流程

public:
	virtual void sample_fitPoints();
	virtual void fit();
	virtual void pia();

protected:
	const int maxIterNum;
	const double eps;
	double error;
	vector<FitPoint2D> fitPoints;
};
#endif // !PIAMETHOD_H

