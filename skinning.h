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
	virtual void parameterize();    // 曲线参数化得到u方向参数 u_knots
	virtual void init() = 0;        // 根据NUUBSCurve初始化T-preimage
	virtual void insert() = 0;      // 在T-preimage中插入中间节点
	virtual void calculate();       // 计算流程

	// update的辅助函数，用于计算基函数的线性组合系数
	void basis_split(const map<double, Node<Point3d>*>& fewer_map, const map<double, Node<Point3d>*>& more_map, map<double, Point3d>& coeff);
	
	// update coordinates of control points by the formula from (nasri 2012)
	// aX' + bW + cY' = V
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
