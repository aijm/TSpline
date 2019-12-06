#include "TsplineVolume.h"


TsplineVolume::TsplineVolume(const TsplineVolume & other)
{
	for (auto entry : other.w_map) {
		Mesh3d* mesh = new Mesh3d(*(entry.second));
		this->w_map[entry.first] = mesh;
	}
	this->w_knots = other.w_knots;
	this->id = -1;
}

TsplineVolume::~TsplineVolume()
{
	for (auto entry : w_map) {
		delete entry.second;
		entry.second = NULL;
	}
}

Point3d TsplineVolume::eval(double u, double v, double w)
{
	Point3d res;

	/*int i = 0;
	for (auto w_it = w_map.begin(); w_it != w_map.end(); w_it++, i++) {
		Mesh3d* tmesh = w_it->second;
		double w_basis = Basis(w_knots, w, i);
		res.add(w_basis * tmesh->eval(u, v));
	}*/
	VectorXd w_knots_change = w_knots;
	w_knots_change(3) = 0.0; w_knots_change(w_knots_change.size() - 4) = 1.0;
	int wid = FindSpan(w_knots_change, w);
	for (int i = wid - 1; i <= wid + 2; i++) {
		double w_basis = Basis(w_knots_change, w, i - 2);
		Mesh3d* mesh = w_map[w_knots(i)];
		res.add(w_basis*mesh->eval(u, v));
	}

	return res;
}

int TsplineVolume::readVolume(string filename)
{
	ifstream in(filename);
	if (!in.is_open()) {
		cout << "failed to open file: " << filename << endl;
		return -1;
	}
	int num = 0;
	in >> num;
	w_knots = Eigen::VectorXd::Zero(num);
	for (int i = 0; i < num; i++) {
		double temp = 0;
		in >> temp;
		w_knots(i) = temp;
	}
	for (int i = 0; i < num - 4; i++) {
		Mesh3d* tmesh = new Mesh3d();
		tmesh->loadMesh(in);
		w_map[w_knots(i + 2)] = tmesh;
	}

	return 0;
}

int TsplineVolume::saveVolume(string filename)
{
	ofstream out(filename + ".vol");
	if (!out.is_open()) {
		cout << "failed to create or open file: " << filename <<".vol"<< endl;
		return -1;
	}
	out << w_knots.size() << endl;
	for (int i = 0; i < w_knots.size(); i++) {
		out << w_knots(i) << " ";
	}
	out << endl;
	for (auto entry : w_map) {
		entry.second->saveMesh(out);
	}
	return 0;
	
}

void TsplineVolume::insert(double w, Mesh3d * mesh)
{
	assert(w > 0.0 && w < 1.0);
	if (w_map.find(w) != w_map.end()) {
		cout << "knots " << w << " in w direction has existed!" << endl;
		exit(1);
	}
	vec_insert(w_knots, w);
	w_map[w] = mesh;
}

int TsplineVolume::get_num() const
{
	int num = 0;
	for (auto entry : w_map) {
		num += entry.second->get_num();
	}
	return num;
}

void TsplineVolume::drawParamCurve()
{
	string filename = "../out/OBJ/edges.obj";
	ofstream out(filename);
	if (!out.is_open()) {
		cout << "failed to open file: " + filename << endl;
		return;
	}

	assert(viewer != NULL);
	const int samples = 5;
	
	Eigen::MatrixXd points;

	for (auto w_it = w_map.begin(); w_it != w_map.end(); w_it++) {
		cout << "----" << endl;
		double w = w_it->first;
		if (w == 0.0 || w == 1.0) {
			continue;
		}
		if (w == 0.0001) {
			w = 0.0;
		}
		if (w == 0.9999) {
			w = 1.0;
		}
		Mesh3d* tmesh = w_it->second;


		auto nodes = tmesh->nodes;

		// nodes of one layer
		for (int i = 0; i < nodes.size(); i++) {
			double u = nodes[i]->s[2];
			double v = nodes[i]->t[2];
			bool valid = (u == 0.0 || u == 1.0) && v != 0.0001 && v != 0.9999;
			valid = valid || ((v == 0.0 || v == 1.0) && u != 0.0001 && u != 0.9999);

			// if a node of current layer has the same (u, v) of the next layer, connect them
			if (valid && next(w_it) != w_map.end()) {
				double w_next = next(w_it)->first; // w value of the next layer
				Mesh3d* tmesh_next = next(w_it)->second; // mesh of the next layer 
				auto node_next = tmesh_next->get_node(u, v);
				if (w_next == 0.9999) {
					w_next = 1.0;
				}
				if (node_next != NULL) {
					vector<Point3d> points(samples + 1);
					out << samples << endl;
					for (int i = 0; i <= samples; i++) {
						double w_val = w + 1.0 * i * (w_next - w) / samples;
						points[i] = eval(u, v, w_val);
						out << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
					}
					for (int i = 0; i < samples; i++) {
						MatrixXd point1 = MatrixXd::Zero(1, 3);
						MatrixXd point2 = MatrixXd::Zero(1, 3);
						point1.row(0) = points[i].toVectorXd().transpose();
						point2.row(0) = points[i + 1].toVectorXd().transpose();
						(*viewer).data().add_edges(point1, point2, black);
					}
					
				}
			}
			

		}
		for (auto iter = tmesh->s_map.begin(); iter != tmesh->s_map.end(); ++iter) {
			if (iter->first == 0.0 || iter->first == 1.0) {
				continue;
			}
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[2]) {
					double s = iter->first;
					if (s == 0.0001) {
						s = 0.0;
					}
					if (s == 0.9999) {
						s = 1.0;
					}
					double t = iter1->first;
					double t_next = (iter1->second)->adj[2]->t[2];
					if (t == 0.0 || t == 0.9999) {
						continue;
					}
					if (t == 0.0001) {
						t = 0.0;
					}
					if (t_next == 0.9999) {
						t_next = 1.0;
					}
					vector<Point3d> points(samples + 1);
					out << samples << endl;
					for (int i = 0; i <= samples; i++) {
						double t_val = t + 1.0 * i * (t_next - t) / samples;
						points[i] = eval(s, t_val, w);
						out << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
					}
					for (int i = 0; i < samples; i++) {
						//cout << "t1, t2: " << t1 << ", " << t2 << endl;
 						MatrixXd point1 = MatrixXd::Zero(1, 3);
						MatrixXd point2 = MatrixXd::Zero(1, 3);
						point1.row(0) = points[i].toVectorXd().transpose();
						point2.row(0) = points[i + 1].toVectorXd().transpose();
						(*viewer).data().add_edges(point1, point2, black);
					}
				}
			}
		}
		for (auto iter = tmesh->t_map.begin(); iter != tmesh->t_map.end(); ++iter) {
			if (iter->first == 0.0 || iter->first == 1.0) {
				continue;
			}
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[1]) {
					double t = iter->first;
					if (t == 0.0001) {
						t = 0.0;
					}
					if (t == 0.9999) {
						t = 1.0;
					}
					double s = iter1->first;
					double s_next = (iter1->second)->adj[1]->s[2];
					if (s == 0.0 || s == 0.9999) {
						continue;
					}
					if (s == 0.0001) {
						s = 0.0;
					}
					if (s_next == 0.9999) {
						s_next = 1.0;
					}
					vector<Point3d> points(samples + 1);
					out << samples << endl;
					for (int i = 0; i <= samples; i++) {
						double s_val = s + 1.0 * i * (s_next - s) / samples;
						points[i] = eval(s_val, t, w);
						out << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
					}
					for (int i = 0; i < samples; i++) {
						double s1 = s + 1.0 * i * (s_next - s) / samples;
						double s2 = s + 1.0 * (i + 1) * (s_next - s) / samples;
						MatrixXd point1 = MatrixXd::Zero(1, 3);
						MatrixXd point2 = MatrixXd::Zero(1, 3);
						point1.row(0) = points[i].toVectorXd().transpose();
						point2.row(0) = points[i + 1].toVectorXd().transpose();
						(*viewer).data().add_edges(point1, point2, black);
					}
					
				}
			}
		}
	}
}

void TsplineVolume::drawTmesh()
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)

	(*viewer).data().add_label(Eigen::Vector3d(-0.02, -0.02, -0.02), "O");
	(*viewer).data().add_label(Eigen::Vector3d(1.05, 0, 0), "u");
	(*viewer).data().add_label(Eigen::Vector3d(0, 1.05, 0), "v");
	(*viewer).data().add_label(Eigen::Vector3d(0, 0, 1.05), "w");

	Eigen::MatrixXd points;

	for (auto w_it = w_map.begin(); w_it != w_map.end(); w_it++) {
		double w = w_it->first;
		Mesh3d* tmesh = w_it->second;
		

		auto nodes = tmesh->nodes;

		// nodes of one layer
		for (int i = 0; i < nodes.size(); i++) {
			Eigen::MatrixXd point(1, 3);
			point << nodes[i]->s[2], nodes[i]->t[2], w;

			(*viewer).data().add_points(point, red);

			points.conservativeResize(points.rows() + 1, 3);
			points.row(points.rows() - 1) = point;

			// 若与下一层有相同的(u,v)节点，则连接成边
			if (next(w_it) != w_map.end()) {
				double w_next = next(w_it)->first; // 下一层的 w 值
				Mesh3d* tmesh_next = next(w_it)->second; // 下一层的网格
				auto node_next = tmesh_next->get_node(point(0, 0), point(0, 1));
				if (node_next != NULL) {
					Eigen::MatrixXd point_next = point;
					point_next(0, 2) = w_next;
					(*viewer).data().add_edges(point, point_next, white);
				}
			}
			
		}
		Eigen::MatrixXd P1 = Eigen::MatrixXd::Zero(1, 3);
		Eigen::MatrixXd P2 = Eigen::MatrixXd::Zero(1, 3);
		P1(0, 2) = w;
		P2(0, 2) = w;
		for (auto iter = tmesh->s_map.begin(); iter != tmesh->s_map.end(); ++iter) {
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[2]) {
					P1(0, 0) = iter->first;
					P1(0, 1) = iter1->first;
					P2(0, 0) = iter->first;
					P2(0, 1) = (iter1->second)->adj[2]->t[2];
					(*viewer).data().add_edges(P1, P2, white);
				}
			}
		}

		for (auto iter = tmesh->t_map.begin(); iter != tmesh->t_map.end(); ++iter) {
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[1]) {
					P1(0, 0) = iter1->first;
					P1(0, 1) = iter->first;
					P2(0, 0) = (iter1->second)->adj[1]->s[2];
					P2(0, 1) = iter->first;
					(*viewer).data().add_edges(P1, P2, white);
				}
			}
		}
	}
	(*viewer).core.align_camera_center(points); // center
}

void TsplineVolume::drawControlpolygon()
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)
	Eigen::MatrixXd points;

	for (auto w_it = w_map.begin(); w_it != w_map.end(); w_it++) {
		double w = w_it->first;
		Mesh3d* tmesh = w_it->second;
		auto nodes = tmesh->nodes;

		// control points of one layer
		for (int i = 0; i < nodes.size(); i++) {
			Eigen::MatrixXd point(1, 3);
			array2matrixd(nodes[i]->data, point);

			(*viewer).data().add_points(point, red);
			points.conservativeResize(points.rows() + 1, 3);
			points.row(points.rows() - 1) = point;

			// 若与下一层有相同的(u,v)节点，则连接成边
			if (next(w_it) != w_map.end()) {
				double w_next = next(w_it)->first; // 下一层的 w 值
				Mesh3d* tmesh_next = next(w_it)->second; // 下一层的网格
				auto node_next = tmesh_next->get_node(nodes[i]->s[2], nodes[i]->t[2]);
				if (node_next != NULL) {
					Eigen::MatrixXd point_next;
					array2matrixd(node_next->data, point_next);
					(*viewer).data().add_edges(point, point_next, yellow);
				}
			}
		}

		Eigen::MatrixXd P1, P2;

		for (auto iter = tmesh->s_map.begin(); iter != tmesh->s_map.end(); ++iter) {
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[2]) {
					array2matrixd(iter1->second->data, P1);
					array2matrixd((iter1->second->adj[2])->data, P2);
					(*viewer).data().add_edges(P1, P2, green);
				}
			}
		}

		for (auto iter = tmesh->t_map.begin(); iter != tmesh->t_map.end(); ++iter) {
			for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if ((iter1->second)->adj[1]) {
					array2matrixd(iter1->second->data, P1);
					array2matrixd((iter1->second->adj[1])->data, P2);
					(*viewer).data().add_edges(P1, P2, blue);
				}
			}
		}
		
	}
	(*viewer).core.align_camera_center(points);
	

}