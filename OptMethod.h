#ifndef OPTMETHOD_H
#define OPTMETHOD_H

#include "PiaMethod.h"

class OptMethod :public Skinning {
public:

	OptMethod(const vector<NURBSCurve>& _curves) :Skinning(_curves) {

	}

	void init() override;		// 根据NUUBSCurve初始化T-preimage
	void insert() override;		// 在T-preimage中插入中间节点
	void calculate() override;  // 计算流程

	static double integral(std::function<double(double, double)> func, double x0 = 0, double x1 = 1, double y0 = 0, double y1 = 1);

private:
	void sample_fitPoints();
	void getM();
	void getN();
	void getB();

	

private:
	vector<FitPoint> fitPoints;
	MatrixXd M; // m_ij
	MatrixXd N; // n_ij
	MatrixXd B; // 

};

#endif // !OPTMETHOD_H

