#include "NurbsPia.h"



void NurbsPia::init()
{
	pia_surface.u_num = num;
	pia_surface.v_num = num;
	pia_surface.u_order = 4;
	pia_surface.v_order = 4;
	pia_surface.dimension = 3;
	pia_surface.isRational = false;
	pia_surface.controlPw = vector<MatrixXd>(num + 1);
	
	VectorXd params(num + 1);
	for (int i = 0; i <= num; i++) {
		pia_surface.controlPw[i] = MatrixXd::Zero(num + 1, 3);
		params[i] = 1.0 * i / num;
	}
	
	VectorXd knots = VectorXd::Zero(num + 5);
	knots(0) = 0.0; knots(1) = 0.0; knots(2) = 0.0; knots(3) = 0.0;
	knots(num + 4) = 1.0; knots(num + 3) = 1.0; knots(num + 2) = 1.0; knots(num + 1) = 1.0;
	knots.segment(4, num - 3) = params.segment(2, num - 3);
	pia_surface.uknots = knots;
	pia_surface.vknots = knots;

	points = vector<MatrixXd>(num + 1, MatrixXd::Zero(num + 1, 3));
	eval = vector<MatrixXd>(num + 1, MatrixXd::Zero(num + 1, 3));
	for (int i = 0; i <= num; i++) {
		for (int j = 0; j <= num; j++) {
			double u = 1.0 * i / num;
			double v = 1.0 * j / num;
			points[j].row(i) = surface.eval(u, v);
		}
	}
	pia_surface.controlPw = points;
}

void NurbsPia::fit()
{
	error = 0.0;
	for (int i = 0; i <= num; i++) {
		for (int j = 0; j <= num; j++) {
			double u = 1.0 * i / num;
			double v = 1.0 * j / num;
			eval[j].row(i) = pia_surface.eval(u, v);
			/*cout << eval[j].row(i) << endl;
			cout << points[j].row(i) << endl;
			cout << "\n\n" << endl;*/
			error += (points[j].row(i) - eval[j].row(i)).norm();
		}
	}
	
	error /= (num + 1);
	error /= (num + 1);
}

void NurbsPia::pia()
{
	for (int i = 0; i < maxIterNum; i++) {
		// 计算差向量并更新曲面控制点
		for (int i = 0; i <= num; i++) {
			for (int j = 0; j <= num; j++) {
				pia_surface.controlPw[j].row(i) += points[j].row(i) - eval[j].row(i);
			}
		}

		fit();
		cout << "iter: " << i + 1 << ", error: " << error << endl;
		if (error < eps) {
			break;
		}

	}
}

NURBSSurface NurbsPia::calculate()
{
	fit();
	pia();
	return pia_surface;
}
