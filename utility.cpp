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


	bool loadpoints(std::string name, Eigen::MatrixXd &mat) {
		ifstream in(name);
		if (!in.is_open()) {
			cout << "error: can't open file: " + name << endl;
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
		//cout << "matrix: \n" << mat << endl;
		return true;
	}
	
	bool savepoints(string name, const Eigen::MatrixXd &mat) {
		name += ".dat";
		ofstream out(name);
		if (!out.is_open()) {
			cout << "error: can't open file: " + name << endl;
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
	void TsplineSimplify(const NURBSSurface & surface, Mesh3d & tspline, int maxIterNum, double eps)
	{
		assert(surface.dimension == 3 && surface.u_order == 4 && surface.v_order == 4);
		// create a tspline format of surface
		Mesh3d origin;
		Eigen::VectorXd u_knots = surface.uknots;
		Eigen::VectorXd v_knots = surface.vknots;
		u_knots(3) = 0.0001; u_knots(u_knots.size() - 4) = 0.9999;
		v_knots(3) = 0.0001; v_knots(v_knots.size() - 4) = 0.9999;
		
		map<double, int> u_map;
		map<double, int> v_map;
		for (int i = 2; i <= u_knots.size()-3; i++) {
			u_map[u_knots(i)] = i;
		}
		for (int i = 2; i <= v_knots.size() - 3; i++) {
			v_map[v_knots(i)] = i;
		}
		
		// 控制顶点对应的参数点的高斯曲率绝对值
		map<double, map<double, double>> curvature;
		double max_curvature = -1;
		double min_curvature = 1e10;

		for (int i = 0; i <= surface.v_num; i++) {
			for (int j = 0; j <= surface.u_num; j++) {
				double u = u_knots(j + 2);
				double v = v_knots(i + 2);
				// 取绝对值后取对数拉伸范围
				double h = log10((abs(surface.guassian_curvature(u, v)) + 1));
				//cout << h << endl;
				curvature[u][v] = h;
				if (h < min_curvature) min_curvature = h;
				if (h > max_curvature) max_curvature = h;
				
				origin.insert_helper(u, v, false);
				auto node = origin.get_node(u, v);
				node->data.fromVectorXd(surface.controlPw[i].row(j));
			}
		}
		cout << "hmin: " << min_curvature << ", " << "hmax : " << max_curvature << endl;
		cout << "origin pool size : " << origin.pool.size() << endl;
		origin.pool.clear();
		if (!origin.check_valid()) {
			cout << "error: invalid tspline mesh!" << endl;
			return;
		}
		//origin.saveMesh("../out/tspline/origin");
		// create an initial tspline patch
		
		Eigen::VectorXd init_knots(4);
		init_knots << 0, 0.0001, 0.9999, 1.0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				double u = init_knots(i);
				double v = init_knots(j);
				tspline.insert_helper(u, v, false);
				tspline.s_map[u][v]->data = origin.s_map[u][v]->data;
			}
		}
		cout << "pool size **********************: " << tspline.pool.size() << endl;
		tspline.pool.clear();
		if (!tspline.check_valid()) {
			cout << "error: invalid tspline mesh!" << endl;
			return;
		}

		// split the rectangle defined by the parameter
		auto split = [&](const std::tuple<double,double,double,double>& rect) {
			// 其实只是分割区域，是拓扑层面的操作，控制点坐标不需要考虑
			int u0 = u_map[get<0>(rect)]; // ulow
			int u1 = u_map[get<1>(rect)]; // uhigh
			int v0 = v_map[get<2>(rect)]; // vlow
			int v1 = v_map[get<3>(rect)]; // vhigh
			
			if (u1 - u0 >= v1 - v0) {
				int u_insert = (u0 + u1) / 2;
				
				tspline.insert_helper(u_knots( u_insert), v_knots(v0),false);
				tspline.merge_all();
				
				tspline.insert_helper(u_knots(u_insert), v_knots(v1),false);
				tspline.merge_all();
			}
			else {
				int v_insert = (v0 + v1) / 2;
				
				tspline.insert_helper(u_knots(u0), v_knots(v_insert),false);
				tspline.merge_all();
				
				tspline.insert_helper(u_knots(u1), v_knots(v_insert),false);
				tspline.merge_all();
			}
		};
		int i = 0;
		while(true) {
			cout << i << " ***********************************************************" << endl;
			cout << "size of nodes: " << tspline.get_num() << endl;
			// 节点插入加细到与B样条曲面一致
			
			Viewer viewer;
			//Mesh3d tspline_copy(tspline);

			//for (auto node : tspline.nodes) {
			//	double u = node->s[2];
			//	double v = node->t[2];
			//	if (origin.get_node(u, v) == 0 || origin.get_node(u,v)->data != node->data) {
			//		cout << "error: *********************" << endl;
			//		return;
			//	}
			//}

			//tspline_copy.setViewer(&viewer);

			//tspline_copy.draw(false, true, true);
			///*for (auto node : origin.nodes) {
			//	MatrixXd P;
			//	array2matrixd(node->data, P);
			//	viewer.data().add_points(P, green);
			//}*/
			//viewer.launch();

			Mesh3d mesh(tspline);
			cout << "size of mesh: " << mesh.get_num() << endl;
			for (auto node : origin.nodes) {
				if (mesh.get_node(node->s[2], node->t[2]) == 0) {
					mesh.insert(node->s[2], node->t[2]);
				}
			}
			
			/*mesh.setViewer(&viewer);
			
			mesh.draw(false, true, true);
			for (auto node : mesh.nodes) {
				auto node_origin = origin.get_node(node->s[2], node->t[2]);
				MatrixXd P1, P2;
				array2matrixd(node->data, P1);
				array2matrixd(node_origin->data, P2);
				viewer.data().add_edges(P1, P2, yellow);
				viewer.data().add_points(P2, green);
			}
			viewer.launch();*/
			
			vector<pair<tuple<double, double, double, double>, double>> regions;
			
			// 先整体计算误差，取出需要split的区域
			for (auto node : mesh.nodes) {

				double u = node->s[2];
				double v = node->t[2];
				double error = (node->data - origin.s_map[u][v]->data).toVectorXd().norm();
				
				double factor = (max_curvature - curvature[u][v]) / (max_curvature - min_curvature);
				if (factor < 0.4) {
					factor = 0.7*factor;
				}
				else {
					factor = (1 + factor) / 2;
				}
				factor = max(factor, 0.05);
		
				if (error < factor*eps) {
					continue;
				}
				auto rects = tspline.region(u, v);
				
				for (auto rect : rects) {
					regions.push_back(make_pair(rect, error - factor*eps));
				}	
			}
			if (regions.empty()) {
				break;
			}
			if (tspline.nodes.size() < origin.nodes.size() * 0.2) {
				// 再按误差排序
				sort(regions.begin(), regions.end(),
					[](const pair<tuple<double, double, double, double>, double>& a,
						const pair<tuple<double, double, double, double>, double>& b) {return a.second > b.second; });
				cout << "最大误差: " << regions[0].second << endl;
				for (auto rect : regions) {
					split(rect.first);
				}
			}
			else {
				// 去除regions的重复元素
				sort(regions.begin(), regions.end(),
					[](const pair<tuple<double, double, double, double>, double>& a,
						const pair<tuple<double, double, double, double>, double>& b) {return a.first < b.first; });
				regions.erase(unique(regions.begin(), regions.end(),
					[](const pair<tuple<double, double, double, double>, double>& a,
						const pair<tuple<double, double, double, double>, double>& b) {return a.first == b.first; }), regions.end());

				// 再按误差排序
				sort(regions.begin(), regions.end(),
					[](const pair<tuple<double, double, double, double>, double>& a,
						const pair<tuple<double, double, double, double>, double>& b) {return a.second > b.second; });

				cout << "total region size: "<<regions.size() << endl;
				// 每次选择误差最大的5个split
				for (int i = 0; i < 10; i++) {
					if (i == regions.size()) {
						break;
					}
					else {
						cout << "误差:" << regions[i].second << endl;
						cout << "split: " << get<0>(regions[i].first)<<", " << get<1>(regions[i].first) << ", "
							<< get<2>(regions[i].first) << ", " << get<3>(regions[i].first) << endl;
						split(regions[i].first);
					}
				}
			}
			


			
			tspline.improve();
			// 更新 tspline

			for (auto node : tspline.nodes) {
				node->data = origin.s_map[node->s[2]][node->t[2]]->data;
			}



			cout << "tspline pool size after split: "<<tspline.pool.size()<<"************************** "<< endl;
			tspline.pool.clear();
			if (!tspline.check_valid()) {
				cout << "error: tspline is not valid!!!" << endl;
			}
			i++;
		}

		cout << "size of nodes: " << tspline.get_num() << endl;
		

	}
};