#include "VolumeSkinning.h"

void VolumeSkinning::parameterize()
{
	/*Eigen::MatrixXd points(surfaces_num, 3);
	for (int i = 0; i < surfaces_num; i++) {
		points.row(i) = centerOfmesh(surfaces[i]).toVectorXd();
	}
	w_params = NURBSCurve::parameterize(points);*/
	w_params = VectorXd::Zero(surfaces_num);
	for (int i = 0; i < surfaces_num; i++) {
		w_params(i) = 1.0*i / (surfaces_num - 1);
	}
	w_params(0) = 0.0001; w_params(w_params.size() - 1) = 0.9999;
	cout << "w_params: " << w_params.transpose() << endl;
}

void VolumeSkinning::init()
{
	volume.w_knots = VectorXd(surfaces_num + 6);
	volume.w_knots(0) = 0; volume.w_knots(1) = 0; volume.w_knots(2) = 0;
	volume.w_knots(surfaces_num+3) = 1; volume.w_knots(surfaces_num + 4) = 1; volume.w_knots(surfaces_num + 5) = 1;
	volume.w_knots.block(3, 0, surfaces_num, 1) = w_params;
	cout << "w_knots: " << volume.w_knots.transpose() << endl;

	volume.w_map[0] = new Mesh3d(surfaces[0]);
	for (int i = 0; i < surfaces_num; i++) {
		volume.w_map[w_params(i)] = new Mesh3d(surfaces[i]);
	}
	volume.w_map[1] = new Mesh3d(surfaces.back());
}

void VolumeSkinning::insert()
{
	assert(surfaces_num >= 2);
	for (int i = 0; i <= surfaces_num - 1; i++) {
		double w_now = w_params(i);
		if (i == 0) {
			double w_insert = 2.0 / 3 * w_now + 1.0 / 3 * w_params(i + 1);
			Mesh3d* mesh = new Mesh3d(surfaces[i]);
			volume.insert(w_insert, mesh);
			for (auto node : mesh->nodes) {
				double u = node->s[2];
				double v = node->t[2];
				node->data = 2.0 / 3 * volume.w_map[w_now]->eval(u, v) + 1.0 / 3 * volume.w_map[w_params(1)]->eval(u, v);
			}
		}
		else if (i == surfaces_num - 1) {
			double w_insert = 2.0 / 3 * w_now + 1.0 / 3 * w_params(i - 1);
			Mesh3d* mesh = new Mesh3d(surfaces[i]);
			volume.insert(w_insert, mesh);
			for (auto node : mesh->nodes) {
				double u = node->s[2];
				double v = node->t[2];
				node->data = 2.0 / 3 * volume.w_map[w_now]->eval(u, v) + 1.0 / 3 * volume.w_map[w_params(i-1)]->eval(u, v);
			}
		}
		else {
			double w_left = w_params(i - 1);
			double w_right = w_params(i + 1);
			double w_insert = 2.0 / 3 * w_now + 1.0 / 3 * w_left;
			Mesh3d* mesh = new Mesh3d(surfaces[i]);
			volume.insert(w_insert, mesh);
			for (auto node : mesh->nodes) {
				double u = node->s[2];
				double v = node->t[2];
				node->data = 2.0 / 3 * volume.w_map[w_now]->eval(u, v) + 1.0 / 3 * volume.w_map[w_left]->eval(u, v);
			}

			w_insert = 2.0 / 3 * w_now + 1.0 / 3 * w_right;
			mesh = new Mesh3d(surfaces[i]);
			volume.insert(w_insert, mesh);
			for (auto node : mesh->nodes) {
				double u = node->s[2];
				double v = node->t[2];
				node->data = 2.0 / 3 * volume.w_map[w_now]->eval(u, v) + 1.0 / 3 * volume.w_map[w_right]->eval(u, v);
			}

		}

	}
}

void VolumeSkinning::calculate()
{
	// 1. compute w_params for surfaces
	parameterize();
	// 2. construct basis 3-dimension T-mesh 
	init();
	// 3. insert intermediate surfaces
	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
	insert();

	// 4. update coordinates of control points by the formula from (nasri 2012)
	// aX' + bW + cY' = V
	update();
}

void VolumeSkinning::update()
{
	map<double, int> id;
	for (int i = 0; i < volume.w_knots.size()-2; i++) {
		id[volume.w_knots(i + 2)] = i;
	}
	for (int i = 1; i < surfaces_num - 1; i++) {
		double w_now = w_params(i);
		int now_id = id[w_now];
		double a = Basis(volume.w_knots, w_now, now_id - 1);
		double b = Basis(volume.w_knots, w_now, now_id);
		double c = Basis(volume.w_knots, w_now, now_id + 1);
		Mesh3d* mesh = volume.w_map[w_now];
		Mesh3d* mesh_left = volume.w_map[volume.w_knots(now_id + 1)];
		Mesh3d* mesh_right = volume.w_map[volume.w_knots(now_id + 3)];

		// aX + bW + cY = V;
		for (int k = 0; k < mesh->nodes.size(); k++) {
			Point3d X = mesh_left->nodes[k]->data;
			Point3d Y = mesh_right->nodes[k]->data;
			Point3d V = surfaces[i].nodes[k]->data;
			mesh->nodes[k]->data = V.add(-a*X).add(-c*Y).scale(1.0 / b);
		}
	}
}

Point3d VolumeSkinning::centerOfmesh(const Mesh3d & mesh)
{
	Point3d res;
	for (auto node : mesh.nodes) {
		res.add(node->data);
	}
	res.scale(1.0 / mesh.nodes.size());
	return res;
}
