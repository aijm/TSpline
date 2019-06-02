#include "Skinning.h"

void Skinning::parameterize()
{
	// 1. compute s-knot for curves
	assert(curves.size() > 0);
	const int dimension = curves[0].controlPw.cols();

	MatrixXd u_cpts(curves_num, dimension);
	for (int i = 0; i < u_cpts.rows(); i++) {
		u_cpts.row(i) = curves[i].controlPw.row(0);
	}
	s_knots = NURBSCurve::parameterize(u_cpts);
	s_knots(0) = 0.0001; s_knots(s_knots.size() - 1) = 0.9999;

	cout << "s_knots: " << s_knots.transpose() << endl;
}

void Skinning::update()
{
	// update coordinates of control points by the formula from (nasri 2012)
	// aX' + bW + cY' = V
	map<double, Point3d> coeff_X; // X'
	map<double, Point3d> coeff_Y; // Y'

	for (int i = 1; i <= curves_num - 2; i++) {
		const double s_now = s_knots(i);
		auto node = tspline.s_map[s_now].begin()->second;
		auto left_node = node->adj[3];
		auto right_node = node->adj[1];
		const double s_left = left_node->s[2];
		const double s_right = right_node->s[2];

		double a = Basis((left_node->s).toVectorXd(), s_now);
		double b = Basis((node->s).toVectorXd(), s_now);
		double c = Basis((right_node->s).toVectorXd(), s_now);

		basis_split(tspline.s_map[s_left], tspline.s_map[s_now], coeff_X); // X'
		basis_split(tspline.s_map[s_right], tspline.s_map[s_now], coeff_Y); // Y'

		for (int j = 0; j <= curves[i].n; j++) {
			double vknot = curves[i].knots(j + 2);
			if (j == 1) {
				vknot = 0.0001;
			}
			if (j == curves[i].n - 1) {
				vknot = 0.9999;
			}
			auto temp_X = coeff_X[vknot];
			auto temp_Y = coeff_Y[vknot];
			temp_X.scale(-a); // -aX'
			temp_Y.scale(-c); // -cY'
			Point3d V;
			V.fromVectorXd(curves[i].controlPw.row(j).transpose());
			tspline.s_map[s_now][vknot]->data = V.add(temp_X).add(temp_Y).scale(1.0 / b); // W=(V-aX'-cY')/b
		}

	}
}

void Skinning::calculate()
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

void Skinning::basis_split(const map<double, Node<Point3d>*>& fewer_map, const map<double, Node<Point3d>*>& more_map, map<double, Point3d>& coeff)
{
	vector<Node<Point3d>> split_pool;
	// 1. copy node of fewer_map to new map
	//vector<Node<T>> fewer_nodes;
	map<double, Node<Point3d>*> new_map;
	for (auto it = fewer_map.begin(); it != fewer_map.end(); ++it) {
		Node<Point3d>* node = new Node<Point3d>();
		*node = *(it->second);
		new_map[it->first] = node;
	}
	for (auto it = more_map.begin(); it != more_map.end(); ++it) {
		if (new_map.find(it->first) == new_map.end()) {
			// split basis function by knot (it->first)
			for (auto it1 = new_map.begin(); it1 != new_map.end(); ++it1) {
				Node<Point3d> temp;
				if (it1->second->split(0, it->first, &temp, true)) {
					split_pool.push_back(temp);
				}
			}
			// merge same node from split_pool to new_map
			for (int i = 0; i < split_pool.size(); i++) {
				double t = split_pool[i].t[2];
				if (new_map.find(t) != new_map.end()) {
					(new_map[t]->data).add(split_pool[i].data);
				}
				else {
					Node<Point3d>* node = new Node<Point3d>();
					*node = split_pool[i];
					new_map[t] = node;
				}
			}
			//bug: split_pool should be cleared
			split_pool.clear();
			cout << "split_pool size: " << split_pool.size() << endl;
		}
	}

	// return data after spliting by reference and delete node
	coeff.clear();
	for (auto it = new_map.begin(); it != new_map.end(); ++it) {
		coeff[it->first] = it->second->data;
		delete it->second;
	}
}
