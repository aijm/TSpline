#include "MinJaeMethod.h"

void MinJaeMethod::init()
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

void MinJaeMethod::insert()
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

void MinJaeMethod::calculate()
{
	// 1. compute s-knot for curves
	parameterize();
	// 2. construct basis T-mesh 
	init();
	// 3. insert intermediate vertices
	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
	insert();

	// 4. 初始化中间点坐标
	inter_init();
	update();
	
	// 5. 迭代更新中间点坐标
	double error = 1.0;
	for (int i = 0; i < maxIterNum; i++) {
		error = inter_update();
		update();
		cout << "*******************iter " << i+1 << ": ,error " << error << "********************" << endl;
		if (error < eps) {
			break;
		}
	}
}

/**
	1. calculate sample points by linear interpolate
	2. fit sample points by B-spline using PIA method
	3. the control points of B-Spline is the initial X,Y
*/
void MinJaeMethod::inter_init()
{
	const int dimension = curves[0].controlPw.cols();
	assert(curves_num >= 2);
	for (int i = 0; i <= curves_num - 2; i++) {
		double s_now = s_knots(i);
		double s_next = s_knots(i + 1);
		auto node_now = (tspline.s_map[s_now].begin())->second;
		auto node_next = (tspline.s_map[s_next].begin())->second;
		double s_inter1 = node_now->adj[1]->s[2];
		double s_inter2 = node_next->adj[3]->s[2];

		MatrixXd sample_inter1(sampleNum + 1, dimension);
		MatrixXd sample_inter2(sampleNum + 1, dimension);

		VectorXd params(sampleNum + 1);

		// calculate sample points by linear interpolate
		for (int j = 0; j <= sampleNum; j++) {
			params(j) = 1.0*j / sampleNum;
			RowVectorXd now_coor = curves[i].eval(params(j));
			RowVectorXd next_coor = curves[i + 1].eval(params(j));
			sample_inter1.row(j) = 2.0 / 3 * now_coor + 1.0 / 3 * next_coor;
			sample_inter2.row(j) = 2.0 / 3 * next_coor + 1.0 / 3 * now_coor;
			
		}
		//(*viewer).data().add_points(sample_inter1, green);
		//(*viewer).data().add_points(sample_inter2, green);

	
		// fit sample points by B-spline using LSPIA with appointed knot vector
		// the control points of B-Spline is the initial X,Y
		VectorXd knots1 = curves[i].knots;
		knots1(3) = 0.0001; knots1(curves[i].n + 1) = 0.9999;

		VectorXd knots2 = curves[i + 1].knots;
		knots2(3) = 0.0001; knots2(curves[i+1].n + 1) = 0.9999;

		MatrixXd cpts1(curves[i].n + 1, 3);
		MatrixXd cpts2(curves[i+1].n + 1, 3);
		for (int j = 0; j <= curves[i].n; j++) {

			cpts1.row(j) = tspline.s_map[s_inter1][knots1(j+2)]->data.toVectorXd();
		}
		for (int j = 0; j <= curves[i+1].n; j++) {
			cpts2.row(j) = tspline.s_map[s_inter2][knots2(j+2)]->data.toVectorXd();
		}
		NURBSCurve inter1, inter2;
		inter1.lspiafit(sample_inter1, params, cpts1, curves[i].knots, 100);
	
		inter2.lspiafit(sample_inter2, params, cpts2, curves[i + 1].knots, 100);

	
		/*inter1.draw(*viewer, false, true, 0.001);
		inter2.draw(*viewer, false, true, 0.001);*/
		
		
		

		inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
		inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;
	
		// update coordinate of intermediate control points and save as inital_cpts
		for (int j = 0; j <= inter1.n; j++) {
			Point3d temp;
	
	
			temp.fromVectorXd(inter1.controlPw.row(j));
	
			initial_cpts[s_inter1][inter1.knots(j + 2)] = temp;
			tspline.s_map[s_inter1][inter1.knots(j + 2)]->data = temp;
	
		}
		for (int j = 0; j <= inter2.n; j++) {
			Point3d temp;
			temp.fromVectorXd(inter2.controlPw.row(j));
	
			initial_cpts[s_inter2][inter2.knots(j + 2)] = temp;
			tspline.s_map[s_inter2][inter2.knots(j + 2)]->data = temp;
		}
	}
	cout << "finished inter_init()!" << endl;
	
}


/**
1. calculate sample points on T-spline surface
2. fit sample points by B-spline using PIA method
3. update coordinate of intermediate control points
*/
double MinJaeMethod::inter_update()
{

	const int dimension = curves[0].controlPw.cols();
	assert(curves_num >= 2);
	map<double, map<double, Point3d>> T_cpts;
	
	for (int i = 0; i <= curves_num - 2; i++) {
		double s_now = s_knots(i);
		double s_next = s_knots(i + 1);
		auto node_now = (tspline.s_map[s_now].begin())->second;
		auto node_next = (tspline.s_map[s_next].begin())->second;
		double s_inter1 = node_now->adj[1]->s[2];
		double s_inter2 = node_next->adj[3]->s[2];
		MatrixXd T_inter1(sampleNum + 1, dimension);
		MatrixXd T_inter2(sampleNum + 1, dimension);

		VectorXd params(sampleNum + 1);

		// 1. calculate sample points on T - spline surface
		for (int j = 0; j <= sampleNum; j++) {
			params(j) = 1.0*j / sampleNum;
			if (j == 0) {
				T_inter1.row(j) = tspline.eval(s_inter1, 0.0001).toVectorXd();
				T_inter2.row(j) = tspline.eval(s_inter2, 0.0001).toVectorXd();
			}
			else if (j == sampleNum) {
				T_inter1.row(j) = tspline.eval(s_inter1, 0.9999).toVectorXd();
				T_inter2.row(j) = tspline.eval(s_inter2, 0.9999).toVectorXd();
			}
			else {
				T_inter1.row(j) = tspline.eval(s_inter1, params(j)).toVectorXd();
				T_inter2.row(j) = tspline.eval(s_inter2, params(j)).toVectorXd();
			}

		}
	
	
		// fit sample points by B-spline using LSPIA with appointed knot vector

		VectorXd knots1 = curves[i].knots;
		knots1(3) = 0.0001; knots1(curves[i].n + 1) = 0.9999;

		VectorXd knots2 = curves[i + 1].knots;
		knots2(3) = 0.0001; knots2(curves[i + 1].n + 1) = 0.9999;

		MatrixXd cpts1(curves[i].n + 1, 3);
		MatrixXd cpts2(curves[i + 1].n + 1, 3);
		for (int j = 0; j <= curves[i].n; j++) {

			cpts1.row(j) = tspline.s_map[s_inter1][knots1(j + 2)]->data.toVectorXd();
		}
		for (int j = 0; j <= curves[i + 1].n; j++) {
			cpts2.row(j) = tspline.s_map[s_inter2][knots2(j + 2)]->data.toVectorXd();
		}

		NURBSCurve inter1, inter2;
		inter1.lspiafit(T_inter1, params, cpts1, curves[i].knots, 100);
	
		inter2.lspiafit(T_inter2, params, cpts2, curves[i + 1].knots, 100);
	
		inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
		inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;
	
		// update coordinate of intermediate control points
		for (int j = 0; j <= inter1.n; j++) {
			T_cpts[s_inter1][inter1.knots(j + 2)].fromVectorXd(inter1.controlPw.row(j).transpose());
		}
	
		for (int j = 0; j <= inter2.n; j++) {
			T_cpts[s_inter2][inter2.knots(j + 2)].fromVectorXd(inter2.controlPw.row(j).transpose());
		}
	}
	
	// update by T_cpts and initial_cpts
	double error = 0.0;
	int count = 0;
	for (auto it = T_cpts.begin(); it != T_cpts.end(); ++it) {
		for (auto it1 = (it->second).begin(); it1 != (it->second).end(); ++it1) {
			count++;
			Point3d delta = (initial_cpts[it->first][it1->first] - (it1->second));
			error += delta.toVectorXd().norm();

			(tspline.s_map[it->first][it1->first]->data).add(delta);
		}
	}
	
	cout << "finished inter_update()!" << endl;
	
	return error / count;
}
