#ifndef MINJAEMETHOD_H
#define MINJAEMETHOD_H

#include "Skinning.h"
class MinJaeMethod : public Skinning {
public:

	MinJaeMethod(const vector<NURBSCurve>& _curves)
		:Skinning(_curves) {

	}

	void init() override;
	void insert() override;
	void calculate() override;

	void inter_update();
};


#endif // !MINJAEMETHOD_H