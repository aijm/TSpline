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
	savepoints("../out/M.mat", M);
	savepoints("../out/N.mat", N);
	savepoints("../out/B.mat", B);

	Eigen::MatrixXd P = (0.00001*M+N).colPivHouseholderQr().solve(B);
	savepoints("../out/P.mat", P);
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
			FitPoint2D point;
			point.param[0] = s_knots(i);
			point.param[1] = 1.0*j / sampleNum;
			point.origin.fromVectorXd(curves[i].eval(point.param[1]).transpose());
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


	for (int i = 0; i <= curves_num - 2; i++) {
		double s_now = s_knots(i);
		double s_next = s_knots(i + 1);
		auto node_now = (tspline.s_map[s_now].begin())->second;
		auto node_next = (tspline.s_map[s_next].begin())->second;
		double s_inter1 = node_now->adj[1]->s[2];
		double s_inter2 = node_next->adj[3]->s[2];

		// calculate sample points by linear interpolate
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint2D point1, point2;
			point1.param[1] = 1.0*j / sampleNum;
			point2.param[1] = point1.param[1];
			point1.param[0] = s_inter1;
			point2.param[0] = s_inter2;

			RowVectorXd now_coor = curves[i].eval(point1.param[1]);
			RowVectorXd next_coor = curves[i + 1].eval(point1.param[1]);
			point1.origin.fromVectorXd(2.0 / 3 * now_coor + 1.0 / 3 * next_coor);
			point2.origin.fromVectorXd(2.0 / 3 * next_coor + 1.0 / 3 * now_coor);
			if (point1.param[1] == 0.0) {
				point1.param[1] = 0.0001;
				point2.param[1] = 0.0001;
			}
			else if (point1.param[1] == 1.0) {
				point1.param[1] = 0.9999;
				point2.param[1] = 0.9999;
			}
			fitPoints.push_back(point1);
			fitPoints.push_back(point2);
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
	assert(M.isApprox(M.transpose()), 1e-5);
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
				double bi = node_i->basis(point.param[0], point.param[1]);
				double bj = node_j->basis(point.param[0], point.param[1]);
				if (i == 83 && bi != 0.0 && bj != 0.0) {
					/*cout << "s: ";
					node_i->s.output(cout);
					cout << endl;

					cout << "t: ";
					node_i->t.output(cout);
					cout << endl;

					cout << "s: ";
					node_j->s.output(cout);
					cout << endl;

					cout << "t: ";
					node_j->t.output(cout);
					cout << endl;

					cout << "u: " << point.param[0] << ", v: " << point.param[1] << endl;
					cout << "bi: " << bi << ", bj: " << bj << endl;
					cout << "\n\n\n" << endl;*/

				}	
				N(i, j) += bi*bj;
			}	
		}
	assert(N.isApprox(N.transpose()), 1e-5);
}

void OptMethod::getB()
{
	const int num = tspline.get_num();
	B = MatrixXd::Zero(num, 3);
	for (int i = 0; i < num; i++) {
		auto node = tspline.get_node(i+1);
		for (const auto& point : fitPoints) {
			Point3d temp = node->basis(point.param[0], point.param[1]) * point.origin;
			for (int j = 0; j < 3; j++) {
				B(i, j) += temp[j];
			}
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
