#ifndef NASRIMETHOD_H
#define NASRIMETHOD_H

#include "Skinning.h"
class NasriMethod : public Skinning {
public:

	NasriMethod(const vector<NURBSCurve>& _curves)
		:Skinning(_curves) {

	}

	void init() override;
	void insert() override;
	void calculate() override;
};


#endif // !NASRIMETHOD_H

