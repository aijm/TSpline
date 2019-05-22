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

}

void MinJaeMethod::inter_update()
{

}
