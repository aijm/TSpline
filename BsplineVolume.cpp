#include "BsplineVolume.h"

Point3d BsplineVolume::eval(double u, double v, double w)
{
	Point3d res;
	// 只需计算 4*4*4个点
	int uid = FindSpan(knot_vector[0], u);
	int vid = FindSpan(knot_vector[1], v);
	int wid = FindSpan(knot_vector[2], w);

	for (int i = uid - 1; i <= uid + 2; i++) {
		for (int j = vid - 1; j <= vid + 2; j++) {
			for (int k = wid - 1; k <= wid + 2; k++) {
				double blend = Basis(knot_vector[0], u, i - 2)*Basis(knot_vector[1], v, j - 2)*Basis(knot_vector[2], w, k - 2);
				res.add(control_grid[i-2][j-2][k-2] * blend);
			}
		}
	}
	return res;
}

int BsplineVolume::readVolume(string filename)
{
	ifstream in(filename);
	if (!in.is_open()) {
		cout << "failed to open file: " << filename << endl;
		return -1;
	}
	
	string sline;
	getline(in, sline);
	getline(in, sline);
	istringstream sin1(sline);
	int x_points, y_points, z_points;
	sin1 >> x_points >> y_points >> z_points;
	control_grid.resize(x_points);
	for (int i = 0; i < x_points; i++)
	{
		control_grid[i].resize(y_points);
		for (int j = 0; j < y_points; j++)
		{
			control_grid[i][j].resize(z_points);
		}
	}
	knot_vector.resize(3);
	knot_vector[0].resize(x_points + 4);
	knot_vector[1].resize(y_points + 4);
	knot_vector[2].resize(z_points + 4);

	getline(in, sline);
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{
				getline(in, sline);
				istringstream sin2(sline);
				sin2 >> control_grid[i][j][k][0] >> control_grid[i][j][k][1] >> control_grid[i][j][k][2];
			}
		}
	}
	getline(in, sline);
	getline(in, sline);
	istringstream sin3(sline);
	for (int i = 0; i < x_points + 4; i++)
	{
		sin3 >> knot_vector[0][i];
	}
	getline(in, sline);
	getline(in, sline);
	istringstream sin4(sline);
	for (int i = 0; i < y_points + 4; i++)
	{
		sin4 >> knot_vector[1][i];
	}
	getline(in, sline);
	getline(in, sline);
	istringstream sin5(sline);
	for (int i = 0; i < z_points + 4; i++)
	{
		sin5 >> knot_vector[2][i];
	}
	in.close();
	cout << filename + " load successful!" << endl;
	return 0;
	
}

int BsplineVolume::saveVolume(string filename)
{
	ofstream out(filename + ".vol");
	if (!out.is_open()) {
		cout << "failed to create or open file: " << filename << ".vol" << endl;
		return -1;
	}
	const int x_num = control_grid.size();
	const int y_num = control_grid[0].size();
	const int z_num = control_grid[0][0].size();

	out << "#resolution of the control grid of BsplineVolume" << endl;
	out << x_num << " " << y_num << " " << z_num << endl;
	out << "#control points of BsplineVolume" << endl;
	for (int i = 0; i < x_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < z_num; k++) {
				out << control_grid[i][j][k][0] << " " << control_grid[i][j][k][1] << " " << control_grid[i][j][k][2] << endl;
			}
		}
	}
	out << "#knot vector in u - direction of BsplineVolume" << endl;
	for (int i = 0; i < knot_vector[0].size(); i++) {
		out << knot_vector[0][i] << " ";
	}
	out << endl;
	out << "#knot vector in v - direction of BsplineVolume" << endl;
	for (int i = 0; i < knot_vector[1].size(); i++) {
		out << knot_vector[1][i] << " ";
	}
	out << endl;
	out << "#knot vector in w - direction of BsplineVolume" << endl;
	for (int i = 0; i < knot_vector[2].size(); i++) {
		out << knot_vector[2][i] << " ";
	}
	out << endl;
	out.close();
	cout << filename + " write successful!" << endl;
	return 0;
}

void BsplineVolume::drawTmesh()
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)

	(*viewer).data().add_label(Eigen::Vector3d(-0.02, -0.02, -0.02), "O");
	(*viewer).data().add_label(Eigen::Vector3d(1.05, 0, 0), "u");
	(*viewer).data().add_label(Eigen::Vector3d(0, 1.05, 0), "v");
	(*viewer).data().add_label(Eigen::Vector3d(0, 0, 1.05), "w");

	

	const int x_num = control_grid.size();
	const int y_num = control_grid[0].size();
	const int z_num = control_grid[0][0].size();

	Eigen::MatrixXd points(x_num*y_num*z_num, 3);
	int count = 0;
	for (int i = 0; i < x_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < z_num; k++){
				
				Eigen::MatrixXd P1(1, 3);
				Eigen::MatrixXd P2(1, 3);
				P1 << knot_vector[0](i + 2), knot_vector[1](j + 2), knot_vector[2](k + 2);
				
				(*viewer).data().add_points(P1, red);
				points.row(count++) = P1;

				if (i != x_num - 1) {
					P2 << knot_vector[0](i + 3), knot_vector[1](j + 2), knot_vector[2](k + 2);
					(*viewer).data().add_edges(P1, P2, white);
				}
				if (j != y_num - 1) {
					knot_vector[0](i + 2), knot_vector[1](j + 3), knot_vector[2](k + 2);
					(*viewer).data().add_edges(P1, P2, white);
				}
				if (k != z_num - 1) {
					knot_vector[0](i + 2), knot_vector[1](j + 2), knot_vector[2](k + 3);
					(*viewer).data().add_edges(P1, P2, white);
				}	
			}
		}
	}
	
	
	(*viewer).core.align_camera_center(points); // center
}

void BsplineVolume::drawControlpolygon()
{
	assert(viewer != NULL); // use setViewer(Viewer* viewer)

	const int x_num = control_grid.size();
	const int y_num = control_grid[0].size();
	const int z_num = control_grid[0][0].size();

	Eigen::MatrixXd points(x_num*y_num*z_num, 3);
	
	int count = 0;
	for (int i = 0; i < x_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < z_num; k++) {
				//cout << i << ", " << j << ", " << k << endl;
				Eigen::MatrixXd P1;
				Eigen::MatrixXd P2;

				array2matrixd(control_grid[i][j][k], P1);
				(*viewer).data().add_points(P1, red);
				points.row(count++) = P1;

				if (i != x_num - 1) {
					array2matrixd(control_grid[i + 1][j][k], P2);
					(*viewer).data().add_edges(P1, P2, white);
				}
				if (j != y_num - 1) {
					array2matrixd(control_grid[i][j + 1][k], P2);
					(*viewer).data().add_edges(P1, P2, white);
				}
				if (k != z_num - 1) {
					array2matrixd(control_grid[i][j][k + 1], P2);
					(*viewer).data().add_edges(P1, P2, white);
				}
			}
		}
	}



	(*viewer).core.align_camera_center(points);

}
