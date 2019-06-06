#pragma once
#include "utility.h"

template<int num>
class FitPoint
{
public:
	double geterror() {
		return (origin - eval).toVectorXd().norm();
	}
	bool inRectangle(const t_mesh::Array<double, num>& low, const t_mesh::Array<double, num>& high) {
		bool res = true;
		for (int i = 0; i < num; i++) {
			if (param[i] < low[i] || param[i] > high[i]) {
				res = false;
			}
		}
		return res;
	}

public:
	Point3d origin; // 要拟合点的坐标
	t_mesh::Array<double, num> param; // 对应参数
	Point3d eval;   // 曲面上对应参数的值 T(u,v) 
	double error;

};

typedef FitPoint<2> FitPoint2D;
typedef FitPoint<3> FitPoint3D;
