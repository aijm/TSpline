#ifndef DRAW_H
#define DRAW_H
#include <igl/opengl/glfw/Viewer.h>

namespace t_mesh {
	const Eigen::RowVector3d   red(1, 0, 0);
	const Eigen::RowVector3d green(0, 1, 0);
	const Eigen::RowVector3d  blue(0, 0, 1);
	const Eigen::RowVector3d white(1, 1, 1);
	const Eigen::RowVector3d black(0, 0, 0);
	const Eigen::RowVector3d deeppink(255, 20, 147);

	

	template<class T,int num>
	void array2matrixd(const Array<T, num> &a, Eigen::MatrixXd &m) {
		m.resize(1, num);
		for (int i = 0; i < num; i++) {
			m(0, i) = 1.0*a[i];
		}
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
		while (i < temp.size() && temp(i)<=t) {
			vec(i) = temp(i);
			i++;
		}
		vec(i) = t;
		while (i < temp.size()) {
			vec(i + 1) = temp(i);
			i++;
		}

	}
	

};

#endif // !DRAW_H

