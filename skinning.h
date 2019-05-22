#ifndef SKINNING_H
#define SKINNING_H

#include "utility.h"
#include "NURBSCurve.h"
using namespace std;
using namespace t_mesh;
class Skinning {
public:
	Skinning(const vector<NURBSCurve>& _curves) :curves(_curves) {
		curves_num = curves.size();
	}
	virtual void parameterize();
	virtual void init() = 0;
	virtual void insert() = 0;
	virtual void calculate() = 0;

	void basis_split(const map<double, Node<Point3d>*>& fewer_map, const map<double, Node<Point3d>*>& more_map, map<double, Point3d>& coeff);
	void update();

	void setViewer(Viewer* _viewer) {
		viewer = _viewer;
	}

public:
	Mesh3d tspline;
protected:
	int curves_num;
	VectorXd s_knots;
	vector<NURBSCurve> curves;
	Viewer* viewer;

};
#endif // !SKINNING_H




//void Skinning_insert1(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
//{
//	const int curves_num = curves.size();

//}
//
//
//
//void intermediate_init(
//	const VectorXd &s_knots,
//	const vector<NURBSCurve>& curves,
//	map<double, map<double, T>>& initial_cpts,
//	Viewer &viewer,
//	int n_sample)
//{
//	const int curves_num = curves.size();
//	const int dimension = curves[0].controlPw.cols();
//	assert(curves_num >= 2);
//	for (int i = 0; i <= curves_num - 2; i++) {
//		double s_now = s_knots(i);
//		double s_next = s_knots(i + 1);
//		auto node_now = (s_map[s_now].begin())->second;
//		auto node_next = (s_map[s_next].begin())->second;
//		double s_inter1 = node_now->adj[1]->s[2];
//		double s_inter2 = node_next->adj[3]->s[2];
//		MatrixXd sample_inter1(n_sample + 1, dimension);
//		MatrixXd sample_inter2(n_sample + 1, dimension);
//		// sample points and linear interpolate
//		for (int j = 0; j <= n_sample; j++) {
//			RowVectorXd now_coor = curves[i].eval(1.0*j / n_sample);
//			RowVectorXd next_coor = curves[i + 1].eval(1.0*j / n_sample);
//			sample_inter1.row(j) = 2.0 / 3 * now_coor + 1.0 / 3 * next_coor;
//			sample_inter2.row(j) = 2.0 / 3 * next_coor + 1.0 / 3 * now_coor;
//		}
//		/*if (true) {
//		viewer.data().add_points(sample_inter1, blue);
//		viewer.data().add_points(sample_inter2, blue);
//		}*/
//
//		/*if (i == 2) {
//		savepoints("sample_inter1" + std::to_string(i) + ".p", sample_inter1);
//		savepoints("sample_inter2" + std::to_string(i) + ".p", sample_inter2);
//		}*/
//
//		// fit the sample points by LSPIA method with appointed knot vector
//		NURBSCurve inter1, inter2;
//		inter1.lspiafit(sample_inter1, curves[i].n + 1, curves[i].knots, 1000);
//
//		inter2.lspiafit(sample_inter2, curves[i + 1].n + 1, curves[i + 1].knots, 1000);
//		/*if (true) {
//		inter1.draw(viewer, true, true, 0.001);
//		inter2.draw(viewer, true, true, 0.001);
//		cout << "inter1.knots: " << inter1.knots.transpose() << endl;
//		cout << "inter2.knots: " << inter2.knots.transpose() << endl;
//		}*/
//
//
//		inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
//		inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;
//
//		// update coordinate of intermediate control points and save as inital_cpts
//		for (int j = 0; j <= inter1.n; j++) {
//			T temp;
//
//
//			temp.fromVectorXd(inter1.controlPw.row(j));
//
//			initial_cpts[s_inter1][inter1.knots(j + 2)] = temp;
//			s_map[s_inter1][inter1.knots(j + 2)]->data = temp;
//
//		}
//		for (int j = 0; j <= inter2.n; j++) {
//			T temp;
//			temp.fromVectorXd(inter2.controlPw.row(j));
//
//			initial_cpts[s_inter2][inter2.knots(j + 2)] = temp;
//			s_map[s_inter2][inter2.knots(j + 2)]->data = temp;
//		}
//	}
//	cout << "finished intermediate_init()!" << endl;
//
//
//}
//
//double intermediate_update(
//	const VectorXd & s_knots,
//	const vector<NURBSCurve>& curves,
//	map<double, map<double, T>>& initial_cpts,
//	Viewer & viewer, int n_sample, bool showiteration)
//{
//	if (showiteration) {
//		draw(viewer, false, false, true);
//	}
//
//	const int curves_num = curves.size();
//	const int dimension = curves[0].controlPw.cols();
//	double error = 0.0;
//	assert(curves_num >= 2);
//	map<double, map<double, T>> T_cpts;
//
//	for (int i = 0; i <= curves_num - 2; i++) {
//		double s_now = s_knots(i);
//		double s_next = s_knots(i + 1);
//		auto node_now = (s_map[s_now].begin())->second;
//		auto node_next = (s_map[s_next].begin())->second;
//		double s_inter1 = node_now->adj[1]->s[2];
//		double s_inter2 = node_next->adj[3]->s[2];
//		MatrixXd T_inter1(n_sample + 1, dimension);
//		MatrixXd T_inter2(n_sample + 1, dimension);
//		// points on T-spline surface
//		for (int j = 0; j <= n_sample; j++) {
//			if (j == 0) {
//				T_inter1.row(j) = (node_now->adj[1]->data).toVectorXd();
//				T_inter2.row(j) = (node_next->adj[3]->data).toVectorXd();
//
//			}
//			else if (j == n_sample) {
//				T_inter1.row(j) = ((s_map[s_now].rbegin())->second->adj[1]->data).toVectorXd();
//				T_inter2.row(j) = ((s_map[s_next].rbegin())->second->adj[3]->data).toVectorXd();
//			}
//			else {
//				T_inter1.row(j) = eval(s_inter1, 1.0*j / n_sample).toVectorXd();
//				T_inter2.row(j) = eval(s_inter2, 1.0*j / n_sample).toVectorXd();
//			}
//
//		}
//		if (showiteration) {
//			viewer.data().add_points(T_inter1, green);
//			viewer.data().add_points(T_inter2, green);
//		}
//
//
//		// fit the sample points by LSPIA method with appointed knot vector
//		NURBSCurve inter1, inter2;
//		inter1.lspiafit(T_inter1, curves[i].n + 1, curves[i].knots, 1000);
//
//		inter2.lspiafit(T_inter2, curves[i + 1].n + 1, curves[i + 1].knots, 1000);
//		if (showiteration) {
//			inter1.draw(viewer, true, true);
//			inter2.draw(viewer, true, true);
//		}
//		/*cout << "inter1.knots: " << inter1.knots.transpose() << endl;
//		cout << "inter2.knots: " << inter2.knots.transpose() << endl;*/
//
//		inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
//		inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;
//
//		// update coordinate of intermediate control points
//		for (int j = 0; j <= inter1.n; j++) {
//			T_cpts[s_inter1][inter1.knots(j + 2)].fromVectorXd(inter1.controlPw.row(j).transpose());
//		}
//
//		for (int j = 0; j <= inter2.n; j++) {
//			T_cpts[s_inter2][inter2.knots(j + 2)].fromVectorXd(inter2.controlPw.row(j).transpose());
//		}
//	}
//
//	// update by T_cpts and initial_cpts
//	for (auto it = T_cpts.begin(); it != T_cpts.end(); ++it) {
//		for (auto it1 = (it->second).begin(); it1 != (it->second).end(); ++it1) {
//			VectorXd delta = (initial_cpts[it->first][it1->first]).toVectorXd() - (it1->second).toVectorXd();
//			if (delta.norm() > error) {
//				error = delta.norm();
//			}
//			T temp;
//			temp.fromVectorXd(delta);
//			(s_map[it->first][it1->first]->data).add(temp);
//		}
//	}
//
//	cout << "finished intermediate_update()!" << endl;
//
//	return error;
//}
//
//template<class T>
//inline void Skinning_intermediate(
//	const VectorXd & s_knots,
//	const vector<NURBSCurve>& curves,
//	Viewer &viewer)
//{
//	// 1. calculate sample points by linear interpolate
//
//	// 2. fit sample points by B-spline using PIA method
//	//    (the control points of B-Spline is the initial X,Y;
//	//     and the B-Spline is the shape to be preserved.)
//	// 3. update V to W by X,Y using formula aX+bW+cY=V
//	// 4. calculate s-curve on T-spline surface, forcing the s-curve to approximate the shape by iteration
//
//	const int curves_num = curves.size();
//	const int dimension = curves[0].controlPw.cols();
//	map<double, map<double, T>> initial_cpts;
//
//	intermediate_init(s_knots, curves, initial_cpts, viewer);
//	Skinning_update_cross(s_knots, curves);
//
//	const int max_iter_num = 0;
//	const double eps = 1e-5;
//	int iter_num = 0;
//	double error = 1.0;
//	while (error > eps&&iter_num < max_iter_num) {
//
//		if (iter_num == max_iter_num - 1) {
//			error = intermediate_update(s_knots, curves, initial_cpts, viewer, 100, false);
//			iter_num++;
//			cout << "*******************iter " << iter_num << ": ,error " << error << "********************" << endl;
//			break;
//		}
//		else {
//			error = intermediate_update(s_knots, curves, initial_cpts, viewer);
//		}
//		//error = intermediate_update(s_knots, curves, initial_cpts, viewer);
//		Skinning_update_cross(s_knots, curves);
//		iter_num++;
//		cout << "*******************iter " << iter_num << ": ,error " << error << "********************" << endl;
//		//Skinning_update_cross(s_knots, curves);
//
//	}
//}
//
//template<class T>
//inline void Skinning_update_cross(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
//{
//	
//}
//// nasri 2012 local T-spline Skinning
//// knot vector of every curve should be [0,0,0,0,...,1,1,1,1]
//// but I need to change multiple knots to single knots [0,0,0,0.0001,...,0.9999,1,1,1]
//template<class T>
//inline void Skinning(const vector<NURBSCurve>& curves, Viewer & viewer)
//{
//	assert(curves.size() > 0);
//	const int curves_num = curves.size();
//	//const int dimension = curves[0].controlPw.cols();
//
//	// 1. compute s-knot for curves
//	VectorXd s_knots = Skinning_parameterize(curves);
//
//	// 2. construct basis T-mesh 
//	Skinning_init(s_knots, curves);
//
//	// 3. insert intermediate vertices
//	// the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
//	//Skinning_insert(s_knots, curves);
//	Skinning_insert1(s_knots, curves);
//	Skinning_intermediate(s_knots, curves, viewer);
//
//	// 4. update coordinates of control points by the formula from (nasri 2012)
//	// aX' + bW + cY' = V
//	//Skinning_update_cross(s_knots, curves);
//}
//
//
//
//template<class T>
//inline void pia_fit(const VectorXd & s_params, const vector<VectorXd>& ttparams, const vector<NURBSCurve>& curves, Viewer &viewer, int max_iterations, double eps)
//{
//	// 每个网格点对应: Q --> 要逼近的点， params --> 对应的参数值
//	map<double, map<double, T>> Q;
//	map<double, map<double, pair<double, double>>> params;
//
//	// 初始化T样条的T-preimage,并将控制点坐标设置为要拟合的点坐标Q
//	VectorXd s_nodes(s_params.size());
//	s_nodes(0) = 0.0; s_nodes(1) = 0.0001;
//	s_nodes(s_nodes.size() - 2) = 0.9999; s_nodes(s_nodes.size() - 1) = 1.0;
//	s_nodes.segment(2, s_nodes.size() - 4) = s_params.segment(2, s_nodes.size() - 4);
//	//cout << "s_nodes:" << s_nodes << endl;
//	for (int i = 0; i < s_nodes.size(); i++) {
//		double s = s_nodes(i);
//		double mapped_s = s_params(i);
//		// 取出曲线节点向量的参数部分
//		/*VectorXd t_params = curves[i].knots.segment(3, curves[i].n - 1);
//		t_params[0] = 0.0001; t_params[t_params.size() - 1] = 0.9999;*/
//		VectorXd t_params = ttparams[i];
//		// 计算相应的T样条节点
//		VectorXd t_nodes(t_params.size());
//		t_nodes(0) = 0.0; t_nodes(1) = 0.0001;
//		t_nodes(t_nodes.size() - 2) = 0.9999; t_nodes(t_nodes.size() - 1) = 1.0;
//		t_nodes.segment(2, t_nodes.size() - 4) = t_params.segment(2, t_nodes.size() - 4);
//
//		for (int j = 0; j < t_params.size(); j++) {
//			double t = t_nodes(j);
//			params[s][t] = make_pair(mapped_s, t_params(j));
//			//cout << "params[" << s << ", " << t << "] = " << mapped_s << ", " << t_params(j) << endl;
//			T point;
//			point.fromVectorXd(curves[i].eval(t_params(j)).row(0).transpose());
//			Q[s][t] = point;
//			MatrixXd P;
//			array2matrixd(point, P);
//			viewer.data().add_points(P, green);
//			insert_helper(s, t, false);
//			Node<T>* node = get_node(s, t);
//			(node->data) = point;
//		}
//	}
//
//	double error = 100.0;
//	int iter_num = 0;
//	// 开始迭代
//	while (error > eps && iter_num < max_iterations) {
//		error = 0.0;
//		for (auto s_it = s_map.begin(); s_it != s_map.end(); s_it++) {
//			for (auto t_it = (s_it->second).begin(); t_it != (s_it->second).end(); t_it++) {
//				double s = s_it->first;
//				double t = t_it->first;
//				T cur = eval(params[s][t].first, params[s][t].second);
//				T delta = Q[s][t] - cur;
//				/*MatrixXd P1, P2;
//				array2matrixd(Q[s][t], P1);
//				array2matrixd(cur, P2);
//
//				viewer.data().add_edges(P1, P2, white);*/
//
//				(s_map[s][t]->data).add(delta);
//				error += delta.toVectorXd().norm();
//			}
//		}
//
//		error /= get_num();
//		iter_num++;
//		//cout << "iter: " << iter_num << ", error: " << error << endl;
//	}
//
//}
//
//template<class T>
//inline double maxoffset(double s, int n, const NURBSCurve& curve, double& max_off)
//{
//	int index = 1;
//	max_off = 0.0;
//	for (int i = 1; i < n; i++) {
//		double t = 1.0*i / n;
//		double offset = (curve.eval(t).row(0).transpose() - eval(s, t).toVectorXd()).norm();
//		if (offset > max_off) {
//			max_off = offset;
//			index = i;
//		}
//	}
//	cout << "maxoffset: " << max_off << endl;
//	return 1.0 * index / n;
//}
//
//void t_mesh::pia_Skinning(const vector<NURBSCurve>& curves, Viewer &viewer, int max_iterations, double eps)
//{
//	assert(curves.size() > 0);
//	const int curves_num = curves.size();
//	//const int dimension = curves[0].controlPw.cols();
//
//	// 1. compute s-knot for curves
//	// 比如 要拟合点的参数: 0, 0.4, 0.6, 0.8, 1.0
//	// 则对应的样条节点 0, 0,           0, 0, 0.4, 0.6, 0.8, 1.0, 1.0,          1.0, 1.0
//	VectorXd s_params = Skinning_parameterize(curves);
//	s_params[0] = 0.0001; s_params[s_params.size() - 1] = 0.9999;
//
//
//	vector<VectorXd> ttparams(curves_num);
//	for (int i = 0; i < ttparams.size(); i++) {
//		// 取出曲线节点向量的参数部分
//		VectorXd t_params = curves[i].knots.segment(3, curves[i].n - 1);
//		t_params[0] = 0.0001; t_params[t_params.size() - 1] = 0.9999;
//		ttparams[i] = t_params;
//	}
//
//	/*for (int i = 0; i < ttparams.size(); i++) {
//	const int n = 10;
//	VectorXd t_params(n + 1);
//	for (int i = 0; i <= n; i++) {
//	t_params(i) = 1.0*i / n;
//	}
//	t_params(0) = 0.0001;
//	t_params(n) = 0.9999;
//	ttparams[i] = t_params;
//	}*/
//
//	for (int i = 0; i < 10; i++) {
//		this->clear();
//		cout << i << ": *****************************" << endl;
//		/*for (int j = 0; j < ttparams.size(); j++) {
//		cout << "curve " << j << ":    " << ttparams[j].transpose() << endl;
//		}*/
//		pia_fit(s_params, ttparams, curves, viewer, max_iterations, eps);
//		cout << "num of nodes: " << get_num() << endl;
//		const int n = 50; // [0,1]区间分为n份
//		for (int i = 0; i < curves_num; i++) {
//			double s = s_params(i);
//			double max_off = 0.0;
//			double t = maxoffset(s, n, curves[i], max_off);
//			//t_mesh::insert(ttparams[i], t);
//
//			if (max_off > eps) {
//				t_mesh::insert(ttparams[i], t);
//			}
//		}
//	}
//}
//
//template<class T>
//inline void basis_split(const map<double, Node<T>*>& fewer_map, const map<double, Node<T>*>& more_map, map<double, T>& coeff)
//{
//	
//
//}
//
