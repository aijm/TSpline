#include "NasriMethod.h"

void NasriMethod::init()
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

void NasriMethod::insert()
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

void NasriMethod::calculate()
{
	// 1. compute s-knot for curves
	parameterize();
	// 2. construct basis T-mesh 
	init();
	// 3. insert intermediate vertices
	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
	insert();

	// 4. update coordinates of control points by the formula from (nasri 2012)
	// aX' + bW + cY' = V
	update();

}
