#include "utility.h"
namespace t_mesh {
	// Blending function N[s0,s1,s2,s3,s4](p)
	double Basis(const Eigen::MatrixXd &knotvector, double t, int i, int p)
	{
		//cout << "knotvector:\n" << knotvector.transpose() << endl;
		//int p = knotvector.size() - 1;
		assert(p >= 1);
		if (p == 1) {
			if (t >= knotvector(i) && t < knotvector(i + 1)) {
				return 1.0;
			}
			else {
				return 0.0;
			}
		}

		double a = knotvector(i + p - 1) - knotvector(i);
		double b = knotvector(i + p) - knotvector(i + 1);
		a = (a == 0.0) ? 0.0 : (t - knotvector(i)) / a;
		b = (b == 0.0) ? 0.0 : (knotvector(i + p) - t) / b;
		return a*Basis(knotvector, t, i, p - 1) + b*Basis(knotvector, t, i + 1, p - 1);
	}


		bool loadpoints(string name, Eigen::MatrixXd &mat) {
			ifstream in(name);
			if (!in) {
				cout << "error: can't open file" + name << endl;
				return false;
			}
			int rows = 0;
			int cols = 0;
			in >> rows >> cols;
			mat = Eigen::MatrixXd(rows, cols);
			for (int i = 0; i < mat.rows(); i++) {
				for (int j = 0; j < mat.cols(); j++) {
					in >> mat(i, j);
				}
			}
			cout << "matrix: \n" << mat << endl;
			return true;
		}
	
		bool savepoints(string name, const Eigen::MatrixXd &mat) {
			ofstream out(name);
			if (!out) {
				cout << "error: can't open file" + name << endl;
				return false;
			}
			out << mat.rows() << " " << mat.cols() << endl;
			out << mat;
			//cout << "matrix: \n" << mat << endl;
			return true;
		}
	
		void insert(Eigen::VectorXd &vec, double t) {
			assert(t > vec(0) && t < vec(vec.size() - 1));
			Eigen::VectorXd temp = vec;
			vec.resize(temp.size() + 1);
			int i = 0;
			while (i < temp.size() && temp(i) <= t) {
				vec(i) = temp(i);
				i++;
			}
			vec(i) = t;
			while (i < temp.size()) {
				vec(i + 1) = temp(i);
				i++;
			}
	
		}
}