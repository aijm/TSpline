#include "Pia.h"

void Pia::init()
{
	Eigen::VectorXd params(4);
	params << 0, 0.0001, 0.9999, 1;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			tspline.insert_helper(params(i), params(j), false);
		}
	}
}

void Pia::adapative_insert()
{

}

void Pia::calculate()
{
	fit();
	pia();
}

void Pia::fit()
{
	error = 0.0;
	for (FitPoint2D& point : fitPoints) {
		point.eval = tspline.eval(point.param[0], point.param[1]);
		point.error = point.geterror();
		error += point.error;
	}
	error /= fitPoints.size();
}

void Pia::pia()
{
	for (int i = 0; i < maxIterNum; i++) {
		// 计算差向量并更新曲面控制点
		for (auto node : tspline.nodes) {
			if (node->s[2] <= 0.0000 || node->s[2] >= 1.0000) {
				continue;
			}
			double sum1 = 0;
			Point3d sum2;
			for (FitPoint2D point : fitPoints) {
				double blend = node->basis(point.param[0], point.param[1]);
				sum1 += blend;
				Point3d delta = point.origin - point.eval;
				delta.scale(blend);
				sum2.add(delta);
			}
			double factor = 0.0;
			if (abs(sum1) > 0.0001) {
				factor = 1.0 / sum1;
			}
			sum2.scale(factor); // 差向量
			node->data.add(sum2); // 更新坐标	
		}

		fit();
		cout << "iter: " << i + 1 << ", error: " << error << endl;
		if (error < eps) {
			break;
		}

	}
}
