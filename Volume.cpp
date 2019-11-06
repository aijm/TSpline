#include "Volume.h"

void Volume::draw(bool tmesh, bool polygon, bool surface, double resolution)
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)
	if (id == -1) {
		id = (*viewer).append_mesh();
	}
	(*viewer).selected_data_index = id;
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

int Volume::saveAsHex(string filename, double resolution)
{
	ofstream out(filename + ".hex");
	if (!out.is_open()) {
		cout << "failed to open file: " + filename + ".hex" << endl;
		return -1;
	}
	ofstream out1(filename + ".jac");
	if (!out1.is_open()) {
		cout << "failed to open file: " + filename + ".jac" << endl;
		return -1;
	}
	// 参数域为 [0,1]*[0,1]*[0,1]
	const int n = 1.0 / resolution;
	const int n1 = n + 1;
	Eigen::MatrixXd V = Eigen::MatrixXd(n1*n1*n1, 3);

	Eigen::MatrixXi Hex = Eigen::MatrixXi(n*n*n, 8);
	
	// jacobian 值
	vector<double> jacobian(Hex.rows() * 8);

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

				Vector3d du = V.row(v_id + 1) - V.row(v_id);
				Vector3d dv = V.row(v_id + n1) - V.row(v_id);
				Vector3d dw = V.row(v_up_id) - V.row(v_id);
				jacobian[hex_id * 8] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = V.row(v_id + n1 + 1) - V.row(v_id + n1);
				dv = -dv;
				dw = V.row(v_up_id + n1) - V.row(v_id + n1);
				jacobian[hex_id * 8 + 1] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = -du;
				dv = V.row(v_id + 1) - V.row(v_id + n1 + 1);
				dw = V.row(v_up_id + n1 + 1) - V.row(v_id + n1 + 1);
				jacobian[hex_id * 8 + 2] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = V.row(v_id) - V.row(v_id + 1);
				dv = -dv;
				dw = V.row(v_up_id + 1) - V.row(v_id + 1);
				jacobian[hex_id * 8 + 3] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id + 1) - V.row(v_up_id);
				dv = V.row(v_up_id + n1) - V.row(v_up_id);
				dw = V.row(v_id) - V.row(v_up_id);
				jacobian[hex_id * 8 + 4] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id + n1 + 1) - V.row(v_up_id + n1);
				dv = -dv;
				dw = V.row(v_id + n1) - V.row(v_up_id + n1);
				jacobian[hex_id * 8 + 5] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = -du;
				dv = V.row(v_up_id + 1) - V.row(v_up_id + n1 + 1);
				dw = V.row(v_id + n1 + 1) - V.row(v_up_id + n1 + 1);
				jacobian[hex_id * 8 + 6] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id) - V.row(v_up_id + 1);
				dv = -dv;
				dw = V.row(v_id + 1) - V.row(v_up_id + 1);
				jacobian[hex_id * 8 + 7] = (du.cross(dv)).normalized().dot(dw.normalized());

				if (reverse) {
					for (int id = 0; id <= 7; id++) {
						jacobian[hex_id * 8 + id] = -jacobian[hex_id * 8 + id];
					}
				}
				/*for (int id = 0; id <= 7; id++) {
					if (jacobian[hex_id * 8 + id] <= 0) {
						cout << "save jacobian: " << jacobian[hex_id * 8 + id] << endl;
						break;
					}
				}*/


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
	sort(jacobian.begin(), jacobian.end());
	out1 << "# jacobian value" << endl;
	out1 << "min: " << jacobian.front() << endl;
	out1 << "max: " << jacobian.back() << endl;
	cout << "minJacobian: " << jacobian.front() << endl;
	cout << "maxJacobian: " << jacobian.back() << endl;

	for (const auto& v : jacobian) {
		out1 << v << endl;
	}
	out.close();
	out1.close();
	return 0;
}

void Volume::drawVolume(double resolution)
{
	
	// 参数域为 [0,1]*[0,1]*[0,1]
	const int n = 1.0 / resolution;
	const int n1 = n + 1;
	Eigen::MatrixXd V = Eigen::MatrixXd(n1*n1*n1, 3);
	Eigen::MatrixXi F = Eigen::MatrixXi(n*n * 12, 3);

	for (int k = 0; k <= n; k++) {
		for (int j = 0; j <= n; j++) {
			for (int i = 0; i <= n; i++) {
				Point3d point = eval(1.0*i / n, 1.0*j / n, 1.0*k / n);
				V.row(k*n1*n1 + j*n1 + i) = point.toVectorXd();
			}
		}
	}
	// 画出8个点
	Eigen::MatrixXd corner(8, 3);
	corner.row(0) = V.row(0);
	corner.row(1) = V.row(n);
	corner.row(2) = V.row(n*n1);
	corner.row(3) = V.row(n*n1+n);
	corner.row(4) = V.row(n*n1*n1);
	corner.row(5) = V.row(n*n1*n1+n);
	corner.row(6) = V.row(n*n1*n1+n*n1);
	corner.row(7) = V.row(n1*n1*n1-1);

	(*viewer).data().add_label(corner.row(0), "000");
	(*viewer).data().add_label(corner.row(1), "100");
	(*viewer).data().add_label(corner.row(2), "010");
	(*viewer).data().add_label(corner.row(3), "110");
	(*viewer).data().add_label(corner.row(4), "001");
	(*viewer).data().add_label(corner.row(5), "101");
	(*viewer).data().add_label(corner.row(6), "011");
	(*viewer).data().add_label(corner.row(7), "111");

	(*viewer).data().add_points(corner, red);

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
	



	Eigen::MatrixXi Hex = Eigen::MatrixXi(n*n*n, 8);

	// jacobian 值
	vector<double> jacobian(Hex.rows() * 8);
	// 六面体网格
	int hex_id = 0;
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				int v_id = k*n1*n1 + j*n1 + i;
				int v_up_id = v_id + n1*n1;
				Hex.row(hex_id) << v_id, v_id + n1, v_id + n1 + 1, v_id + 1,
					v_up_id, v_up_id + n1, v_up_id + n1 + 1, v_up_id + 1;

				Vector3d du = V.row(v_id + 1) - V.row(v_id);
				Vector3d dv = V.row(v_id + n1) - V.row(v_id);
				Vector3d dw = V.row(v_up_id) - V.row(v_id);
				jacobian[hex_id * 8] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = V.row(v_id + n1 + 1) - V.row(v_id + n1);
				dv = -dv;
				dw = V.row(v_up_id + n1) - V.row(v_id + n1);
				jacobian[hex_id * 8 + 1] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = -du;
				dv = V.row(v_id + 1) - V.row(v_id + n1 + 1);
				dw = V.row(v_up_id + n1 + 1) - V.row(v_id + n1 + 1);
				jacobian[hex_id * 8 + 2] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = V.row(v_id) - V.row(v_id + 1);
				dv = -dv;
				dw = V.row(v_up_id + 1) - V.row(v_id + 1);
				jacobian[hex_id * 8 + 3] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id + 1) - V.row(v_up_id);
				dv = V.row(v_up_id + n1) - V.row(v_up_id);
				dw = V.row(v_id) - V.row(v_up_id);
				jacobian[hex_id * 8 + 4] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id + n1 + 1) - V.row(v_up_id + n1);
				dv = -dv;
				dw = V.row(v_id + n1) - V.row(v_up_id + n1);
				jacobian[hex_id * 8 + 5] = (du.cross(dv)).normalized().dot(dw.normalized());

				du = -du;
				dv = V.row(v_up_id + 1) - V.row(v_up_id + n1 + 1);
				dw = V.row(v_id + n1 + 1) - V.row(v_up_id + n1 + 1);
				jacobian[hex_id * 8 + 6] = (dv.cross(du)).normalized().dot(dw.normalized());

				du = V.row(v_up_id) - V.row(v_up_id + 1);
				dv = -dv;
				dw = V.row(v_id + 1) - V.row(v_up_id + 1);
				jacobian[hex_id * 8 + 7] = (du.cross(dv)).normalized().dot(dw.normalized());


				if (reverse) {
					for (int id = 0; id <= 7; id++) {
						jacobian[hex_id * 8 + id] = -jacobian[hex_id * 8 + id];
					}
				}

				for (int id = 0; id <= 7; id++) {
					if (jacobian[hex_id * 8 + id] <= 0) {
						//cout << "negative jacobian: " << jacobian[hex_id * 8 + id] << endl;
						MatrixXd p0 = MatrixXd::Zero(1, 3);
						p0.row(0) = V.row(Hex(hex_id, 0));
						MatrixXd p1 = MatrixXd::Zero(1, 3);
						p1.row(0) = V.row(Hex(hex_id, 1));
						MatrixXd p2 = MatrixXd::Zero(1, 3);
						p2.row(0) = V.row(Hex(hex_id, 2));
						MatrixXd p3 = MatrixXd::Zero(1, 3);
						p3.row(0) = V.row(Hex(hex_id, 3));
						MatrixXd p4 = MatrixXd::Zero(1, 3);
						p4.row(0) = V.row(Hex(hex_id, 4));
						MatrixXd p5 = MatrixXd::Zero(1, 3);
						p5.row(0) = V.row(Hex(hex_id, 5));
						MatrixXd p6 = MatrixXd::Zero(1, 3);
						p6.row(0) = V.row(Hex(hex_id, 6));
						MatrixXd p7 = MatrixXd::Zero(1, 3);
						p7.row(0) = V.row(Hex(hex_id, 7));
						
						(*viewer).data().add_edges(p0, p1, white);
						(*viewer).data().add_edges(p1, p2, white);
						(*viewer).data().add_edges(p2, p3, white);
						(*viewer).data().add_edges(p0, p3, white);
						(*viewer).data().add_edges(p4, p5, white);
						(*viewer).data().add_edges(p5, p6, white);
						(*viewer).data().add_edges(p6, p7, white);
						(*viewer).data().add_edges(p7, p4, white);
						(*viewer).data().add_edges(p0, p4, white);
						(*viewer).data().add_edges(p1, p5, white);
						(*viewer).data().add_edges(p2, p6, white);
						(*viewer).data().add_edges(p3, p7, white);


						
						/*for (int m = 0; m < 8; m++) {
							MatrixXd point = MatrixXd::Zero(1, 3);
							point.row(0) = V.row(Hex(hex_id, m));
							(*viewer).data().add_points(point, white);
						}*/
						break;
					}
				}

				hex_id++;
			}
		}
	}

	(*viewer).data().set_mesh(V, F);

}
