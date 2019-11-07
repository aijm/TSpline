#include "PiaMinJaeMethod.h"

std::tuple<double,double,double,double> PiaMinJaeMethod::param_helper_points()
{
	const int n = 1.0 / 0.01;
	const int n1 = n + 1;
	Eigen::MatrixXd V = Eigen::MatrixXd(n1*n1, 3);
	Eigen::MatrixXi F = Eigen::MatrixXi(2 * n*n, 3);
	// discretize T-Spline Surface into triangular mesh(V,F) in libigl mesh structure
	// calculate 

	for (int j = 0; j <= n; j++)
		for (int i = 0; i <= n; i++) {
			Eigen::MatrixXd curvePoint;
			double u = 1.0*i / n;
			double v = 1.0*j / n;
			array2matrixd(tspline.eval(u, v), curvePoint);
			V.row(j*n1 + i) = curvePoint;
		}

	for (int j = 0; j<n; j++)
		for (int i = 0; i < n; i++) {
			int V_index = j*(n + 1) + i;
			int F_index = 2 * j*n + 2 * i;
			F.row(F_index) << V_index, V_index + 1, V_index + n + 1;
			F.row(F_index + 1) << V_index + n + 1, V_index + 1, V_index + n + 2;
		}
	// 根据顶点索引返回对应参数
	auto eval_param_uv = [n](Vector3i face){
		double u = 0, v = 0;
		for (int k = 0; k < 3; k++) {
			int ii = face(k) % (n + 1);
			int jj = face(k) / (n + 1);
			u += 1.0*ii / n;
			v += 1.0*jj / n;
		}
		return std::make_tuple(u / 3, v / 3);
	};

	// 计算点到三角网格的投影，据此计算对于参数值
	MatrixXd query_points(helper_points.size(), 3);
	
	for (int i = 0; i < helper_points.size(); i++) {
		query_points.row(i) = helper_points[i].origin.toVectorXd();
	}
	//(*viewer).data().add_points(query_points, red);

	VectorXd sqrD(query_points.rows()); // 点到三角网格的最短距离平方
	VectorXi face(query_points.rows()); // 点到三角网格最短距离点所在面 在F中的索引
	MatrixXd closest_point(query_points.rows(), 3); // 最近的点的坐标

	igl::point_mesh_squared_distance(query_points, V, F, sqrD, face, closest_point);
	// 画出结果
	/*(*viewer).data().add_edges(query_points, closest_point, green);
	for (int i = 0; i < query_points.rows(); i++) {
		MatrixXd v0(1, 3), v1(1, 3), v2(1, 3);
		cout << "0,1,2: " << F(face(i), 0) << ", " << F(face(i), 1) << ", " << F(face(i), 2) << endl;
		v0 = V.row(F(face(i), 0));
		v1 = V.row(F(face(i), 1));
		v2 = V.row(F(face(i), 2));
		(*viewer).data().add_edges(v0, v1, red);
		(*viewer).data().add_edges(v0, v2, red);
		(*viewer).data().add_edges(v1, v2, red);
	}*/

	// 计算对应参数值, 为投影所在三角面片三个顶点的参数值的重心
	double umin = 1.0, vmin = 1.0, umax = 0.0, vmax = 0.0;
	for (int i = 0; i < query_points.rows(); i++) {
		auto param = eval_param_uv(F.row(face(i)));
		helper_points[i].param[0] = std::get<0>(param);
		helper_points[i].param[1] = std::get<1>(param);
		double u = helper_points[i].param[0];
		double v = helper_points[i].param[1];
		if (u < umin) umin = u;
		if (u > umax) umax = u;
		if (v < vmin) vmin = v;
		if (v > vmax) vmax = v;

		cout << "i,u,v: " <<i <<", "<< u << ", " << v << endl;
	}
	return std::make_tuple(1.2*umin, 1.2*umax, 1.2*vmin, 1.2*vmax);
}

void PiaMinJaeMethod::init()
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

void PiaMinJaeMethod::insert()
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

void PiaMinJaeMethod::calculate()
{
	parameterize();
	init();
	insert();
	sample_fitPoints_2();
	fitPoints = curve_points;
	fitPoints.insert(fitPoints.end(), inter_points.begin(), inter_points.end());
	fit();
	pia();
	update();
	
	for (int i = 0; i < 10; i++) {
		pia();
		update();
	}

	for (int j = 0; j < 5; j++) {
		auto grid = param_helper_points(); // 辅助点参数化
		t_mesh::Array<double, 2> low;
		t_mesh::Array<double, 2> high;
		low[0] = std::get<0>(grid); high[0] = std::get<1>(grid);
		low[1] = std::get<2>(grid); high[1] = std::get<3>(grid);
		//fitPoints = curve_points;
		fitPoints.clear();
		fitPoints.insert(fitPoints.end(), helper_points.begin(), helper_points.end());
		for (auto point : inter_points) {
			if (!point.inRectangle(low,high)) {
				fitPoints.push_back(point);
			}
		}
		for (auto point : curve_points) {
			if (!point.inRectangle(low,high)) {
				fitPoints.push_back(point);
			}
		}
		fit();
		pia();
	}
	update();
	
}

void PiaMinJaeMethod::sample_fitPoints_1()
{
	const int sampleNum = 10;
	for (int i = 1; i < curves_num - 1; i++) {
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
			curve_points.push_back(point);

			/*MatrixXd P;
			array2matrixd(point.origin, P);
			(*viewer).data().add_points(P, green);*/
		}
	}

	// 直接从nurbs蒙皮结果上采样
	NURBSSurface skin;
	VectorXd params = s_knots;
	params(0) = 0; params(params.size() - 1) = 1;
	Viewer viewer;
	skin.skinning(curves, params, viewer);
	const int num = 10;
	for (int i = 0; i <= num; i++) {
		bool valid = true;
		for (int k = 0; k < s_knots.size(); k++) {
			if (abs(s_knots(k) - 1.0*i / num) < 0.01) {
				valid = false;
				break;
			}
		}
		if (!valid) {
			continue;
		}
		for (int j = 0; j <= num; j++) {
			FitPoint2D point;
			point.param[0] = 1.0*i / num;
			
			point.param[1] = 1.0*j / num;
			point.origin.fromVectorXd(skin.eval(point.param[0], point.param[1]));
			inter_points.push_back(point);
		}
	}

}

void PiaMinJaeMethod::sample_fitPoints_2()
{
	const int sampleNum = 50;
	for (int i = 1; i < curves_num - 1; i++) {
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
			curve_points.push_back(point);

			/*MatrixXd P;
			array2matrixd(point.origin, P);
			(*viewer).data().add_points(P, red);*/
		}
	}

	// 纵向采样，拟合出一个B样条曲线
	// others --> 10 , 10
	// helicoidal --> 20, 20
	const int v_sample_num = 20;
	const int u_sample_num = 20;

	VectorXd params = s_knots;
	params(0) = 0; params(params.size() - 1) = 1;
	cout << "params: \n" << params << endl;
	VectorXd knots(params.size() + 6);
	knots(0) = 0; knots(1) = 0, knots(2) = 0;
	knots(knots.size()-1) = 1; knots(knots.size() - 2) = 1, knots(knots.size() - 3) = 1;
	knots.block(3, 0, params.size(), 1) = params;
	cout << "knots for interpolate: \n" << knots << endl;

	vector<NURBSCurve> sample_curves(v_sample_num + 1);
	for (int i = 0; i <= v_sample_num; i++) {
		double v = 1.0*i / v_sample_num;
		MatrixXd points(curves_num, 3);
		for (int j = 0; j < curves_num; j++) {
			points.row(j) = curves[j].eval(v);
		}
		sample_curves[i].interpolate(points, knots);
		//sample_curves[i].draw(*viewer, false, true, 0.001, blue);
	}


	for (int i = 0; i <= u_sample_num; i++) {
		bool valid = true;
		for (int k = 0; k < s_knots.size(); k++) {
			if (abs(s_knots(k) - 1.0*i / u_sample_num) < 0.01) {
				valid = false;
				break;
			}
		}
		if (!valid) {
			continue;
		}
		for (int j = 0; j <= v_sample_num; j++) {
			FitPoint2D point;
			point.param[0] = 1.0*i / u_sample_num;
			point.param[1] = 1.0*j / v_sample_num;

			point.origin.fromVectorXd(sample_curves[j].eval(point.param[0]));
			inter_points.push_back(point);

			/*MatrixXd P;
			array2matrixd(point.origin, P);
			(*viewer).data().add_points(P, green);*/
		}
	}
	
}

void PiaMinJaeMethod::sample_fitPoints()
{
	const int sampleNum = 50;
	for (int i = 1; i < curves_num - 1; i++) {
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint2D point;
			point.param[0] = s_knots(i);
			point.param[1] = 1.0*j / sampleNum;
			point.origin.fromVectorXd(curves[i].eval(point.param[1]).transpose());
			/*if (point.param[1] == 0.0) {
				point.param[1] = 0.0001;
			}
			else if (point.param[1] == 1.0) {
				point.param[1] = 0.9999;
			}*/
			curve_points.push_back(point);

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
			/*if (point1.param[1] == 0.0) {
				point1.param[1] = 0.0001;
				point2.param[1] = 0.0001;
			}
			else if (point1.param[1] == 1.0) {
				point1.param[1] = 0.9999;
				point2.param[1] = 0.9999;
			}*/
			inter_points.push_back(point1);
			inter_points.push_back(point2);

			/*MatrixXd P;
			array2matrixd(point1.origin, P);
			(*viewer).data().add_points(P, green);
			array2matrixd(point2.origin, P);
			(*viewer).data().add_points(P, green);*/
		}
	}
	cout << "number of points: " << curve_points.size() + inter_points.size() << endl;
}

void PiaMinJaeMethod::pia()
{
	for (int i = 0; i < maxIterNum; i++) {
		// 计算差向量并更新曲面控制点
		for (auto node : tspline.nodes) {
			if (node->s[2] <= 0.0001 || node->s[2] >= 0.9999) {
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
