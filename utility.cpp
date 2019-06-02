#include "utility.h"
namespace t_mesh {
	// the index of parameter t in knot vector
	// example knots=[0, 0, 0, 0, 0.4, 0.6, 1, 1, 1, 1]; t= 0.5, n = 5, p = 3 
	// return 4
	int FindSpan(const Eigen::MatrixXd &knots, double t, int p) {
		const int n = knots.size() - p - 2;
		if (t == knots(n + 1)) return n;
		int low = p;
		int high = n+1;
		assert(t >= knots(low) && t < knots(high));

		int mid = (low + high) / 2;
		while (t < knots(mid) || t >= knots(mid + 1))
		{
			if (t < knots(mid)) high = mid;
			else low = mid;
			mid = (low + high) / 2;
		}
		return mid;
	}

	// Blending function N[s0,s1,s2,s3,s4](t)
	double Basis(const Eigen::MatrixXd &knots, double t, int i, int p)
	{
		const int m = knots.size() - 1;
		// 特殊情况
		if (i == 0) {
			int count = 0;
			for (int j = 0; j <= p; j++) {
				if (abs(knots(j) - t) <= 0.0001) {
					count++;
				}
			}
			if (count == p + 1) {
				return 1.0;
			}
		}
		if (i == m - p - 1) {
			int count = 0;
			for (int j = 0; j <= p; j++) {
				if (abs(knots(i+j+1) - t) <= 0.0001) {
					count++;
				}
			}
			if (count == p + 1) {
				return 1.0;
			}
		}
		/*if ((i == 0 && t == knots(0)) ||
			(i == m - p - 1 && t == knots(m))) {
			return 1.0;
		}*/
		// 根据局部性
		if (t < knots(i) || t >= knots(i + p + 1)) {
			return 0.0;
		}
		Eigen::VectorXd N(p + 1);
		// 初始化0次的基函数
		for (int j = 0; j <= p; j++) {
			if (t >= knots(i + j) && t < knots(i + j + 1)) N(j) = 1.0;
			else N(j) = 0.0;
		}
		//cout << "N: \n" << N << endl;

		// 计算三角形表
		for (int k = 1; k <= p; k++) {
			//cout << "k=1:" << endl;
			double saved = 0.0;
			if (N(0) == 0.0) saved = 0.0; 
			else saved = (t - knots(i)) * N(0) / (knots(i + k) - knots(i));

			for (int j = 0; j < p - k + 1; j++) {
				double Uleft = knots(i + j + 1);
				double Uright = knots(i + j + k + 1);
				if (N(j + 1) == 0.0) { 
					N(j) = saved;
					saved = 0.0;
				}
				else {
					double temp = N(j + 1) / (Uright - Uleft);
					N(j) = saved + (Uright - t) * temp;
					saved = (t - Uleft) * temp;
				}
				//cout << N(j) << endl;
			}
		}
		return N(0);
	}

	// Berivative of Blending function N[s0,s1,s2,s3,s4](t)
	Eigen::RowVectorXd DersBasis(const Eigen::MatrixXd &knots, double t, int i, int p)
	{
		if (t == 0.0) {
			t = 0.0001;
		}
		if (t == 1.0) {
			t = 0.9999;
		}
		const int m = knots.size() - 1;
		Eigen::RowVectorXd ders = Eigen::RowVectorXd::Zero(p + 1); // k阶导数, k= 0,1,2,...,p
	
		// 根据局部性
		if (t < knots(i) || t >= knots(i + p + 1)) {
			for (int k = 0; k <= p; k++) ders(k) = 0.0;
			return ders;
		}
		
		Eigen::MatrixXd N = Eigen::MatrixXd::Zero(p + 1, p + 1);
		// 初始化0次的基函数
		for (int j = 0; j <= p; j++) {
			if (t >= knots(i + j) && t < knots(i + j + 1)) N(j, 0) = 1.0;
			else N(j, 0) = 0.0;
		}

		// 计算三角形表
		for (int k = 1; k <= p; k++) {
			double saved = 0.0;
			if (N(0,k-1) == 0.0) saved = 0.0;
			else saved = (t - knots(i)) * N(0,k-1) / (knots(i + k) - knots(i));

			for (int j = 0; j < p - k + 1; j++) {
				double Uleft = knots(i + j + 1);
				double Uright = knots(i + j + k + 1);
				if (N(j + 1, k - 1) == 0.0) {
					N(j, k) = saved;
					saved = 0.0;
				}
				else {
					double temp = N(j + 1, k - 1) / (Uright - Uleft);
					N(j, k) = saved + (Uright - t) * temp;
					saved = (t - Uleft) * temp;
				}
			}
		}
		//cout << "N: \n" << N << endl;
		ders(0) = N(0, p); // 函数值

		// 计算导数
		for (int k = 1; k <= p; k++) {
			Eigen::VectorXd ND = N.col(p - k); // 载入正确的列

			for (int jj = 1; jj <= k; jj++) {
				double saved = 0.0;
				if (ND(0) == 0.0) saved = 0.0;
				else saved = ND(0) / (knots(i + p - k + jj) - knots(i));

				for (int j = 0; j < k - jj + 1; j++) {
					double Uleft = knots(i + j + 1);
					double Uright = knots(i + j + p - k + jj + 1);
					if (ND(j + 1) == 0.0) {
						ND(j) = (p - k + jj) * saved;
						saved = 0.0;
					}
					else {
						double temp = ND(j + 1) / (Uright - Uleft);
						ND(j) = (p - k + jj)*(saved - temp);
						saved = temp;
					}
				}
			}
			ders(k) = ND(0);
		}
		return ders;
	}


	// Blending function N[s0,s1,s2,s3,s4](p)
	double Basis1(const Eigen::MatrixXd &knotvector, double t, int i, int p)
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
	
		void vec_insert(Eigen::VectorXd &vec, double t) {
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