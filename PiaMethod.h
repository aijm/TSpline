#ifndef PIAMETHOD_H
#define PIAMETHOD_H

#include "skinning.h"
class PiaMethod : public Skinning {
public:
	PiaMethod(const vector<NURBSCurve>& _curves)
		:Skinning(_curves) {

	}
	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 在T-preimage中插入中间节点
	void calculate() override;  // 计算流程
};
#endif // !PIAMETHOD_H

