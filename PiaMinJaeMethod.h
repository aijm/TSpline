#ifndef PIAMINJAEMETHOD_H
#define PIAMINJAEMETHOD_H

#include "PiaMethod.h"
class PiaMinJaeMethod : public PiaMethod {
public:
	PiaMinJaeMethod(const vector<NURBSCurve>& _curves, int _maxIterNum = 100, double _eps = 1e-5)
		:PiaMethod(_curves, _maxIterNum, _eps) {

	}
	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 按一定规则在误差大的地方插入节点，局部加细
	void calculate() override;  // 计算流程

public:
	void sample_fitPoints() override;
};

#endif // !PIAMINJAEMETHOD_H

