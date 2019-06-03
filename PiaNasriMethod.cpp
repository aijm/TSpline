#include "PiaNasriMethod.h"

void PiaNasriMethod::init()
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

void PiaNasriMethod::insert()
{
	// 3. insert intermediate vertices
	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
	assert(curves_num >= 3);
	for (int i = 0; i <= curves_num - 2; i++) {
		double s_now = s_knots(i);
		auto s_nodes = tspline.s_map[s_now];
		for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
			if (it->second->adj[1]) {
				double s_mid = (s_now + s_knots(i + 1)) / 2;
				tspline.insert_helper(s_mid, it->first, false);
				auto node = tspline.get_node(s_mid, it->first);
				//(node->data).add(it->second->data);
				//(node->data).add(it->second->adj[1]->data);
				// bug, it->second->adj[1] has changed to current node
				(node->data).add(node->adj[1]->data);
				(node->data).add(node->adj[3]->data);
				(node->data).scale(0.5);
			}
		}
	}
	tspline.pool.clear();
	if (!tspline.check_valid()) {
		cout << "skinning: invalid T-mesh!" << endl;
		return;
	}
}

void PiaNasriMethod::calculate()
{
	parameterize();
	init();
	insert();
	sample_fitPoints();
	fit();
	pia();
	update();
}

void PiaNasriMethod::sample_fitPoints()
{
	const int sampleNum = 10;
	for (int i = 1; i < curves_num-1; i++) {
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint point;
			point.u = s_knots(i);
			point.v = 1.0*j / sampleNum;
			point.origin.fromVectorXd(curves[i].eval(point.v).transpose());
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


	for (int i = 0; i <= curves_num - 2; i++) {
		double s_now = s_knots(i);
		double s_next = s_knots(i + 1);
		auto node_now = (tspline.s_map[s_now].begin())->second;
		auto node_next = (tspline.s_map[s_next].begin())->second;
		double s_inter1 = node_now->adj[1]->s[2];
		double s_inter2 = node_next->adj[3]->s[2];

		// calculate sample points by linear interpolate
		for (int j = 0; j <= sampleNum; j++) {
			FitPoint point1, point2;
			point1.v = 1.0*j / sampleNum;
			point2.v = point1.v;
			point1.u = s_inter1;
			point2.u = s_inter2;

			RowVectorXd now_coor = curves[i].eval(point1.v);
			RowVectorXd next_coor = curves[i + 1].eval(point1.v);
			point1.origin.fromVectorXd(2.0 / 3 * now_coor + 1.0 / 3 * next_coor);
			point2.origin.fromVectorXd(2.0 / 3 * next_coor + 1.0 / 3 * now_coor);
			if (point1.v == 0.0) {
				point1.v = 0.0001;
				point2.v = 0.0001;
			}
			else if (point1.v == 1.0) {
				point1.v = 0.9999;
				point2.v = 0.9999;
			}
			fitPoints.push_back(point1);
			fitPoints.push_back(point2);
		}
	}
	cout << "number of points: " << fitPoints.size() << endl;
}

void PiaNasriMethod::pia()
{
	for (int i = 0; i < maxIterNum; i++) {
		// 计算差向量并更新曲面控制点
		for (auto node : tspline.nodes) {
			if (node->s[2] <= 0.0001 || node->s[2] >= 0.9999) {
				continue;
			}
			double sum1 = 0;
			Point3d sum2;
			for (FitPoint point : fitPoints) {
				double blend = node->basis(point.u, point.v);
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
