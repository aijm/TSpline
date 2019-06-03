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

	// 参数域为 [0,1]*[0,1]*[0,1]
	const int n = 1.0 / resolution;
	const int n1 = n + 1;
	Eigen::MatrixXd V = Eigen::MatrixXd(n1*n1*n1, 3);

	Eigen::MatrixXi Hex = Eigen::MatrixXi(n*n*n, 8);

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

	(*viewer).data().set_mesh(V, F);
	

}
