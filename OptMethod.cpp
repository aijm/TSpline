#include "OptMethod.h"

void OptMethod::init()
{
	// 2. construct basis T-mesh 
	//add 0 and 1
	for (int j = 0; j <= curves[0].n; j++) {
		double vknot = curves[0].knots(j + 2);
		if (j == 1) { vknot = 0.0001; }
		if (j == curves[0].n - 1) { vknot = 0.9999; }

		tspline.insert_helper(0.0, vknot, false);
		auto node = tspline.get_node(0.0, vknot);
		(node->data).fromVectorXd(curves[0].controlPw.row(j));
		/*(node->data).output(cout);
		cout << endl;*/
	}

	for (int i = 0; i < curves_num; i++) {
		for (int j = 0; j <= curves[i].n; j++) {
			double vknot = curves[i].knots(j + 2);
			if (j == 1) { vknot = 0.0001; }
			if (j == curves[i].n - 1) { vknot = 0.9999; }
			tspline.insert_helper(s_knots(i), vknot, false);
			auto node = tspline.get_node(s_knots(i), vknot);
			(node->data).fromVectorXd(curves[i].controlPw.row(j));
			/*(node->data).output(cout);
			cout << endl;*/
			//merge_all();
		}
	}
	for (int j = 0; j <= curves[curves_num - 1].n; j++) {
		double vknot = curves[curves_num - 1].knots(j + 2);
		if (j == 1) { vknot = 0.0001; }
		if (j == curves[curves_num - 1].n - 1) { vknot = 0.9999; }
		tspline.insert_helper(1.0, vknot, false);
		auto node = tspline.get_node(1.0, vknot);
		(node->data).fromVectorXd(curves[curves_num - 1].controlPw.row(j));

		/*(node->data).output(cout);
		cout << endl;*/
	}

	cout << "pool size:" << tspline.pool.size() << endl;
	tspline.pool.clear();

	if (!tspline.check_valid()) {
		cout << "skinning: invalid T-mesh!" << endl;
		return;
	}
}

void OptMethod::insert()
{
	// 3. insert intermediate vertices
	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
	assert(curves_num >= 3);
	for (int i = 0; i <= curves_num - 1; i++) {
		double s_now = s_knots(i);
		auto s_nodes = tspline.s_map[s_now];
		if (i == 0) {
			for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
				double s_insert = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i + 1);
				tspline.insert_helper(s_insert, it->first, false);
				MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i + 1].eval(it->first);
				(tspline.s_map[s_insert][it->first]->data).fromVectorXd(position.row(0).transpose());
			}
		}
		else if (i == curves_num - 1) {
			for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
				double s_insert = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i - 1);
				tspline.insert_helper(s_insert, it->first, false);
				MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i - 1].eval(it->first);
				(tspline.s_map[s_insert][it->first]->data).fromVectorXd(position.row(0).transpose());
			}
		}
		else {
			for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
				double s_left = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i - 1);
				double s_right = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i + 1);
				tspline.insert_helper(s_left, it->first, false);
				MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i - 1].eval(it->first);
				(tspline.s_map[s_left][it->first]->data).fromVectorXd(position.row(0).transpose());

				tspline.insert_helper(s_right, it->first, false);
				position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i + 1].eval(it->first);
				(tspline.s_map[s_right][it->first]->data).fromVectorXd(position.row(0).transpose());
			}
		}

	}
	tspline.pool.clear();
	
	if (!tspline.check_valid()) {
		cout << "Skinning: invalid T-mesh!" << endl;
		return;
	}
}

void OptMethod::calculate()
{
	parameterize();
	init();
	insert();
	sample_fitPoints();
	getM();
	getN();
	getB();
	savepoints("../M.mat", M);
	savepoints("../N.mat", N);
	savepoints("../B.mat", B);

	Eigen::MatrixXd P = N.colPivHouseholderQr().solve(B);
	savepoints("../P.mat", P);
	for (int i = 0; i < tspline.get_num(); i++) {
		auto node = tspline.get_node(i + 1);
		node->data.fromVectorXd(P.row(i).transpose());
	}
	

}

void OptMethod::sample_fitPoints()
{
	const int sampleNum = 100;
	for (int i = 0; i < curves_num; i++) {
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint point;
			point.u = s_knots(i);
			point.v = 1.0*j / sampleNum;
			point.origin.fromVectorXd(curves[i].eval(point.v).row(0).transpose());
			if (point.v == 0.0) {
				point.v = 0.0001;
			}
			else if (point.v == 1.0) {
				point.v = 0.9999;
			}
			fitPoints.push_back(point);

			/*MatrixXd P;
			array2matrixd(point.origin, P);
			(*viewer).data().add_points(P, green);*/
		}
	}
	cout << "number of points: " << fitPoints.size() << endl;
}

void OptMethod::getM()
{
	const int num = tspline.get_num();
	M = MatrixXd::Zero(num, num);
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			auto node_i = tspline.get_node(i+1);
			auto node_j = tspline.get_node(j+1);

			auto lambda = [this,node_i,node_j](double u, double v)->double {
				double res = 0;
				double t0 =  tspline.du2(node_i, u, v)*tspline.du2(node_j, u, v);
				double t1 =  2 * tspline.duv(node_i, u, v)*tspline.duv(node_j, u, v);
				double t2 =  tspline.dv2(node_i, u, v)*tspline.dv2(node_j, u, v);
				res += t0 + t1 + t2;
				/*cout << "t0: " << t0 << endl;
				cout << "t1: " << t1 << endl;
				cout << "t2: " << t2 << endl;*/
				return res;
			};
			M(i, j) = integral(lambda);
		}
	}
	assert(M.isApprox(M.transpose()));
}

void OptMethod::getN()
{
	const int num = tspline.get_num();
	N = MatrixXd::Zero(num, num);

	for(int i = 0;i < N.rows();i++)
		for (int j = 0; j < N.cols(); j++) 
		{
			auto node_i = tspline.get_node(i+1);
			auto node_j = tspline.get_node(j+1);
			for (const auto& point : fitPoints) {
				N(i, j) += node_i->basis(point.u, point.v) * node_j->basis(point.u, point.v);
			}	
		}
	assert(N.isApprox(N.transpose()));
}

void OptMethod::getB()
{
	const int num = tspline.get_num();
	B = MatrixXd::Zero(num, 3);
	for (int i = 0; i < num; i++) {
		auto node = tspline.get_node(i+1);
		for (const auto& point : fitPoints) {
			B.row(i) += node->basis(point.u, point.v) * point.origin.toVectorXd().transpose();
		}
	}
}

double OptMethod::integral(std::function<double(double, double)> func, double x0, double x1,double y0,double y1)
{
	// Gaussian quadrature 计算 func(u,v)在[a,b][a,b]上的积分
	const int n = 5;
	VectorXd w(n);
	w << 0.2369268851, 0.4786286705, 0.5688888888, 0.4786286705, 0.2369268851;
	VectorXd x(n);
	x << -0.9061798459, -0.5384693101, 0.0000000000, 0.5384693101, 0.9061798459;

	double sum = 0;
	double px0 = (x0 + x1) / 2;
	double px1 = (x1 - x0) / 2;

	double py0 = (y0 + y1) / 2;
	double py1 = (y1 - y0) / 2;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum += w(i)*w(j)*func(px0 + px1*x(i), py0 + py1*x(j));
		}
	}
	sum *= px1 * py1;
	return sum;
}
