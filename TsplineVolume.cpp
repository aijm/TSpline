#include "TsplineVolume.h"

void TsplineVolume::draw(bool tmesh, bool polygon, bool surface, double resolution)
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)

	if (tmesh) {
		drawTmesh();
		return;
	}
	if (polygon) {
		drawControlpolygon();
	}
	if (surface) {
		drawVolume(resolution);
	}
}

Point3d TsplineVolume::eval(double u, double v, double w)
{
	Point3d res;
	int i = 0;
	for (auto w_it = w_map.begin(); w_it != w_map.end(); w_it++,i++) {
		Mesh3d* tmesh = w_it->second;
		double w_basis = Basis(w_knots, w, i);
		res.add(w_basis * tmesh->eval(u, v));
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
	out << w_knots.transpose() << endl;
	for (auto entry : w_map) {
		entry.second->saveMesh(out);
	}
	return 0;
	
}

int TsplineVolume::saveAsHex(string filename, double resolution)
{
	ofstream out(filename + ".hex");
	if (!out.is_open()) {
		cout << "failed to open file: " + filename + ".hex" << endl;
		return -1;
	}

	// 参数域为 [0,1]*[0,1]*[0,1]
	const int n = 1.0 / resolution;
	const int n1 = n + 1;
	V = Eigen::MatrixXd(n1*n1*n1, 3);

	Hex = Eigen::MatrixXi(n*n*n, 8);

	for (int k = 0; k <= n; k++) {
		for (int j = 0; j <= n; j++) {
			for (int i = 0; i <= n; i++) {
				Point3d point = eval(1.0*i / n, 1.0*j / n, 1.0*k / n);
				V.row(k*n1*n1 + j*n1 + i) = point.toVectorXd();
			}
		}
	}
	
	// 六面体网格
	int hex_id = 0;
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				int v_id = k*n1*n1 + j*n1 + i;
				int v_up_id = v_id + n1*n1;
				Hex.row(hex_id) << v_id, v_id + n1, v_id + n1 + 1, v_id + 1,
					v_up_id, v_up_id + n1, v_up_id + n1 + 1, v_up_id + 1;

				hex_id++;
			}
		}
	}
	out << "# num of vertices: " << V.rows() << endl;
	out << "# num of Hexahedron: " << Hex.rows() << endl;
	for (int i = 0; i < V.rows(); i++) {
		out << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << endl;
	}
	for (int i = 0; i < Hex.rows(); i++) {
		out << "h ";
		for (int j = 0; j < 8; j++) {
			out << Hex(i, j) << " ";
		}
		out << endl;
	}

	return 0;
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

			points.conservativeResize(point.rows() + 1, 3);
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

void TsplineVolume::drawVolume(double resolution)
{
	// 参数域为 [0,1]*[0,1]*[0,1]
	const int n = 1.0 / resolution;
	const int n1 = n + 1;
	V = Eigen::MatrixXd(n1*n1*n1, 3);
	F = Eigen::MatrixXi(n*n * 12, 3);

	for (int k = 0; k <= n; k++) {
		for (int j = 0; j <= n; j++) {
			for (int i = 0; i <= n; i++) {
				Point3d point = eval(1.0*i / n, 1.0*j / n, 1.0*k / n);
				V.row(k*n1*n1 + j*n1 + i) = point.toVectorXd();
			}
		}
	}
	// 画出六个面
	int f_id = 0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			// bottom
			int v_id = j*n1 + i;
			//int f_id = 2 * (2 * j*n + 2 * i);
			F.row(f_id) << v_id, v_id + n1, v_id + 1;
			F.row(f_id + 1) << v_id + n1, v_id + n1 + 1, v_id + 1;
			// top
			f_id += 2;
			v_id += n*n1*n1;

			F.row(f_id) << v_id, v_id + 1, v_id + n1;
			F.row(f_id + 1) << v_id + 1, v_id + n1 + 1, v_id + n1;
			// front
			f_id += 2;
			v_id = j*n1*n1 + i;
			
			F.row(f_id) << v_id, v_id + 1, v_id + n1*n1;
			F.row(f_id + 1) << v_id + 1, v_id + n1*n1 + 1, v_id + n1*n1;
			// back
			f_id += 2;
			v_id += n*(n + 1);
			
			F.row(f_id) << v_id, v_id + n1*n1, v_id + 1;
			F.row(f_id + 1) << v_id + n1*n1, v_id + n1*n1 + 1, v_id + 1;
			// left
			f_id += 2;
			v_id = j*n1*n1 + i*n1;
			
			F.row(f_id) << v_id, v_id + n1*n1, v_id + n1;
			F.row(f_id + 1) << v_id + n1*n1, v_id + n1*n1 + n1, v_id + n1;
			// right
			f_id += 2;
			v_id += n;
			
			F.row(f_id) << v_id, v_id + n1, v_id + n1*n1;
			F.row(f_id + 1) << v_id + n1, v_id + n1*n1 + n1, v_id + n1*n1;
			f_id += 2;
		}
	}

	(*viewer).data().set_mesh(V, F);
}
