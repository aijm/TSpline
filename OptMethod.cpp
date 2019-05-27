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

}
