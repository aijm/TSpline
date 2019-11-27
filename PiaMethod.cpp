#include "PiaMethod.h"

void PiaMethod::init()
{
	// 曲线上采样得到待拟合的目标点
	
	// 构造初始T-preimage
	// 2. construct basis T-mesh 
	//add 0, 0.0001(重节点)
	for (int j = 0; j <= curves[0].n; j++) {
		double vknot = curves[0].knots(j + 2);
		if (j == 1) { vknot = 0.0001; }
		if (j == curves[0].n - 1) { vknot = 0.9999; }

		tspline.insert_helper(0.0, vknot, false);
		auto node = tspline.get_node(0.0, vknot);
		(node->data).fromVectorXd(curves[0].controlPw.row(j));

		tspline.insert_helper(0.0001, vknot, false);
		node = tspline.get_node(0.0001, vknot);
		(node->data).fromVectorXd(curves[0].controlPw.row(j));

		/*(node->data).output(cout);
		cout << endl;*/
	}
	// add 1, 0.9999(重节点)
	for (int j = 0; j <= curves[curves_num - 1].n; j++) {
		double vknot = curves[curves_num - 1].knots(j + 2);
		if (j == 1) { vknot = 0.0001; }
		if (j == curves[curves_num - 1].n - 1) { vknot = 0.9999; }
		tspline.insert_helper(1.0, vknot, false);
		auto node = tspline.get_node(1.0, vknot);
		(node->data).fromVectorXd(curves[curves_num - 1].controlPw.row(j));

		tspline.insert_helper(0.9999, vknot, false);
		node = tspline.get_node(0.9999, vknot);
		(node->data).fromVectorXd(curves[curves_num - 1].controlPw.row(j));
	}

	cout << "pool size:" << tspline.pool.size() << endl;
	tspline.pool.clear();

	if (!tspline.check_valid()) {
		cout << "skinning: invalid T-mesh!" << endl;
		return;
	}
}

void PiaMethod::insert()
{
	int id = 0;
	for (int i = 0; i < fitPoints.size(); i++) {
		if (fitPoints[i].error > fitPoints[id].error) {
			id = i;
		}
	}
	auto node = tspline.get_knot(fitPoints[id].param[0], fitPoints[id].param[1]);
	cout << "maxnode " << node.s[2] << ", " << node.t[2] << endl;
	double u = (node.s[1] + node.s[3]) / 2.0;
	double v = (node.t[1] + node.t[3]) / 2.0;
	tspline.insert(u, v);
	cout << "insert " << u << ", " << v << endl;
}

void PiaMethod::calculate()
{
	parameterize();
	init();
	sample_fitPoints();
	fit();

	for(int i = 0;i < 5;i++) {
		pia();
		cout << "\n\n\n\n\n";
		insert();
		cout << "\n\n\n\n\n";
		fit();
	}
	
}

void PiaMethod::sample_fitPoints()
{
	const int sampleNum = 10;
	for (int i = 1; i < curves_num-1; i++) {
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint2D point;
			point.param[0] = s_knots(i);
			point.param[1] = 1.0*j / sampleNum;
			point.origin.fromVectorXd(curves[i].eval(point.param[1]).row(0).transpose());
			if (point.param[1] == 0.0) {
				point.param[1] = 0.0001;
			}
			else if (point.param[1] == 1.0) {
				point.param[1] = 0.9999;
			}
			fitPoints.push_back(point);

			/*MatrixXd P;
			array2matrixd(point.origin, P);
			(*viewer).data().add_points(P, green);*/
		}
	}
	cout << "number of points: " << fitPoints.size() << endl;
}

// 计算曲面上对应参数处的坐标 T(u, v)
void PiaMethod::fit()
{
	error = 0.0;
	for (FitPoint2D& point : fitPoints) {
		point.eval = tspline.eval(point.param[0], point.param[1]);
		point.error = point.geterror();
		error += point.error * point.error;
	}
	//error /= fitPoints.size();
	//error = sqrt(error);
}

void PiaMethod::pia()
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
		//cout << "iter: " << i + 1 << ", error: " << error << endl;
		cout << error << endl;
		if (error < eps) {
			break;
		}
		
	}
	

}

