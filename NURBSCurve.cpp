// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "NURBSCurve.h"

/*input format:
_n       : P_0,P_1,...,P_n; _n is the final index 
_k       : order of BSpline
_controlP: P_0,P_1,...,P_n; (n+1) by 2 or 3
_knots   : t_0,t_1,...,t_(n+k); */
NURBSCurve::NURBSCurve(int _n, int _k, MatrixXd _controlP, VectorXd _knots, bool _isRational)
{
	assert(_k >= 1 && _n >= _k - 1 && _controlP.rows() == _n + 1 && _knots.size() == _n + _k + 1);

	n = _n;
	k = _k;
	knots = _knots;
	controlPw = _controlP;
	isRational = _isRational;

}

// load
bool NURBSCurve::loadNURBS(string name){
	ifstream in(name.c_str());
	if(!in){
		return false;
	}
	in >> isRational;
	in >> n >> k;
	int dimension=3;
	in >> dimension;

	controlPw = MatrixXd(n+1,dimension);
	knots = VectorXd(n+k+1);
	// to do.
	for(int i=0;i < controlPw.rows();i++){
		for(int j=0;j< controlPw.cols();j++){
			in >> controlPw(i,j);
		}
	}
	for(int i=0;i<knots.size();i++){
		in >> knots(i);
	}
	return true;
}
// save
// nurbs 格式示例
/*
0              isRational
6 4            n k
3              点的维数
2    3  0      控制点坐标
1.5  1  0
- 1  1  0
- 2   0  0
- 1 - 1  0
1 - 1  0
2    0  0
0    0    0 0 0.25    0.5 0.75    1 1    1    1  节点向量
*/
bool NURBSCurve::saveNURBS(string name){
	if(isRational){
		name += ".cptw";
	}else{
		name += ".cpt";
	}
	ofstream out(name.c_str());
	if(!out){
		return false;
	}
	out<< isRational<<endl;
	out<< n <<" "<< k<<endl;
	out<< controlPw.cols()<<endl;
	out<< controlPw <<endl;
	out<< knots.transpose();
	return true;
}

// find the knot interval of t by binary searching
int NURBSCurve::find_ind(double t) const
{
	
	if (t == knots(n + 1)) return n;
	int low = 0;
	int high = n+k;
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

// evaluate the coordinate of curvePoint with parameter t  
MatrixXd NURBSCurve::eval(double t) const
{
	// find the knot interval of t by binary searching
	int L = find_ind(t); //[t_L,t_(L+1)] 
	assert(L >= k - 1 && L <= n + 1);

	// P_(L-k+1),..,P_L control the interval [t_L,t_(L+1)]
	MatrixXd temp = controlPw.block(L - k + 1, 0, k, controlPw.cols());

	// de-Boor algorithm
	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}

	return temp.row(0);
}

VectorXd NURBSCurve::parameterize(const MatrixXd & points)
{
	
	int K = points.rows() - 1;
	// chord length parametrization, u_0,...,u_K
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i - 1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;
	return params;
}



double NURBSCurve::Basis(const VectorXd & _knots, double _t, int _i, int _p)
{
	const int m = _knots.size() - 1;
	// 特殊情况
	if (_i == 0) {
		int count = 0;
		for (int j = 0; j <= _p; j++) {
			if (abs(_knots(j) - _t) <= 0.0001) {
				count++;
			}
		}
		if (count == _p + 1) {
			return 1.0;
		}
	}
	if (_i == m - _p - 1) {
		int count = 0;
		for (int j = 0; j <= _p; j++) {
			if (abs(_knots(_i + j + 1) - _t) <= 0.0001) {
				count++;
			}
		}
		if (count == _p + 1) {
			return 1.0;
		}
	}
	
	// 根据局部性
	if (_t < _knots(_i) || _t >= _knots(_i + _p + 1)) {
		return 0.0;
	}
	Eigen::VectorXd N(_p + 1);
	// 初始化0次的基函数
	for (int j = 0; j <= _p; j++) {
		if (_t >= _knots(_i + j) && _t < _knots(_i + j + 1)) N(j) = 1.0;
		else N(j) = 0.0;
	}
	//cout << "N: \n" << N << endl;

	// 计算三角形表
	for (int k = 1; k <= _p; k++) {
		//cout << "k=1:" << endl;
		double saved = 0.0;
		if (N(0) == 0.0) saved = 0.0;
		else saved = (_t - _knots(_i)) * N(0) / (_knots(_i + k) - _knots(_i));

		for (int j = 0; j < _p - k + 1; j++) {
			double Uleft = _knots(_i + j + 1);
			double Uright = _knots(_i + j + k + 1);
			if (N(j + 1) == 0.0) {
				N(j) = saved;
				saved = 0.0;
			}
			else {
				double temp = N(j + 1) / (Uright - Uleft);
				N(j) = saved + (Uright - _t) * temp;
				saved = (_t - Uleft) * temp;
			}
			//cout << N(j) << endl;
		}
	}
	return N(0);
}

void NURBSCurve::interpolate(const MatrixXd &points)
{
	// points: P_0, ..., P_K
	//assert(points.rows() >= 4 && (points.cols() == 2 || points.cols() == 3));
	int K = points.rows() - 1;
	// chord length parametrization, u_0,...,u_K
	VectorXd params(points.rows());
	params(0) = 0.0;
	for (int i = 1; i <= K; i++) {
		params(i) = params(i-1) + (points.row(i) - points.row(i - 1)).norm();
	}
	params = params / params(K);
	params(K) = 1.0;

	knots = VectorXd::Zero(K + 7); // u_0,u_0,u_0, u_0,...,u_K, u_K,u_K,u_K
	
	knots.block(3, 0, K + 1, 1) = params;
	knots(K + 6) = 1.0; knots(K + 5) = 1.0; knots(K + 4) = 1.0;
	//cout << "knots:\n" << knots << endl;
	isRational = false;
	if (points.cols() == 4) {
		isRational = true;
	}
	n = K + 2; // X_0,X_1,...,X_(K+2)
	k = 4; // order 4
	controlPw = MatrixXd::Zero(K + 3, points.cols());
	controlPw.row(0) = points.row(0);
	controlPw.row(K + 2) = points.row(K);

	
	//// solve linear equation: AX = b
	MatrixXd X(K + 1, points.cols()); //X_1,...,X_(K+1)
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);
	
	MatrixXd b = points;

	//A(0, 0) = A(K, K) = 1.0;
	double delta0 = params(1) - params(0);
	double delta1 = params(2) - params(1);
	double delta2 = params(K) - params(K - 1);
	double delta3 = params(K - 1) - params(K - 2);
	A(0, 0) = (2.0*delta0 + delta1) / (delta0 + delta1);
	A(0, 1) = -delta0 / (delta0 + delta1);
	A(K, K - 1) = -delta2 / (delta2 + delta3);
	A(K, K) = (2.0*delta2 + delta3) / (delta2 + delta3);


	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = Basis(knots, params(i), i);
		A(i, i) = Basis(knots, params(i), i+1);
		A(i, i + 1) = Basis(knots, params(i), i + 2);

	}
	//cout << "A:\n" << A << endl;
	X = A.inverse()*b;
	//cout << "X:\n" << X << endl;
	controlPw.block(1, 0, K + 1, points.cols()) = X;
	

}

// interpolate with appointed knot vector
void NURBSCurve::interpolate(const MatrixXd &points, const VectorXd &knotvector)
{
	// points: P_0, ..., P_K
	//assert(points.rows() >= 4 && (points.cols() == 2 || points.cols() == 3));
	int K = points.rows() - 1;
	assert(knotvector.size() == points.rows() + 6);
	
	knots = knotvector;
	VectorXd params = knots.block(3, 0, K + 1, 1);

	//cout << "knots:\n" << knots << endl;
	isRational = false;
	if (points.cols() == 4) {
		isRational = true;
	}
	n = K + 2; // X_0,X_1,...,X_(K+2)
	k = 4; // order 4
	controlPw = MatrixXd::Zero(K + 3, points.cols());
	controlPw.row(0) = points.row(0);
	controlPw.row(K + 2) = points.row(K);

	// solve linear equation: AX = b
	MatrixXd X(K + 1, points.cols()); //X_1,...,X_(K+1)
	MatrixXd A = MatrixXd::Zero(K + 1, K + 1);

	MatrixXd b = points;
	

	//A(0, 0) = A(K, K) = 1.0;
	double delta0 = params(1) - params(0);
	double delta1 = params(2) - params(1);
	double delta2 = params(K) - params(K - 1);
	double delta3 = params(K - 1) - params(K - 2);
	A(0, 0) = (2.0*delta0 + delta1) / (delta0 + delta1);
	A(0, 1) = -delta0 / (delta0 + delta1);
	A(K, K - 1) = -delta2 / (delta2 + delta3);
	A(K, K) = (2.0*delta2 + delta3) / (delta2 + delta3);

	for (int i = 1; i <= K - 1; i++) {
		A(i, i - 1) = Basis(knots, params(i), i);
		A(i, i) = Basis(knots, params(i), i+1);
		
		A(i, i + 1) = Basis(knots, params(i), i+2);

	}
	//cout << "A:\n" << A << endl;
	X = A.inverse()*b;
	//cout << "X:\n" << X << endl;
	controlPw.block(1, 0, K + 1, points.cols()) = X;
	
}

void NURBSCurve::piafit(const MatrixXd &points, int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->n = points.rows() - 1;
	this->k = 4;
	this->isRational = false;
	VectorXd params(points.rows());
	params = parameterize(points);
	cout << "params: " << params.transpose() << endl;
	knots = VectorXd(n + k + 1);
	knots(0) = 0.0; knots(1) = 0.0; knots(2) = 0.0; knots(3) = 0.0;
	knots(n + 4) = 1.0; knots(n + 3) = 1.0; knots(n + 2) = 1.0; knots(n + 1) = 1.0;
	knots.segment(4, n - 3) = params.segment(2, n - 3);
	cout << "knots:" << knots.transpose() << endl;

	controlPw = points;
	double error = 1.0;
	int iter_num = 0;
	MatrixXd delta(points.rows(), points.cols());
	while (error>eps && iter_num<max_iter_num) {
		
		for (int j = 0; j <= n; j++) {
			delta.row(j) = points.row(j) - eval(params(j));
		}
		
		controlPw += delta;
		
		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
		cout << "iter: " << iter_num << ", error: " << error << endl;
	}


}

void NURBSCurve::piafit(const MatrixXd &points, const VectorXd &knotvector,int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->n = points.rows() - 1;
	this->k = 4;
	this->isRational = false;
	VectorXd params(points.rows());
	params = parameterize(points);
	knots = knotvector;


	controlPw = points;
	double error = 1.0;
	int iter_num = 0;
	MatrixXd delta(points.rows(), points.cols());
	while (error>eps && iter_num<max_iter_num) {

		for (int j = 0; j <= n; j++) {
			delta.row(j) = points.row(j) - eval(params(j));
		}

		controlPw += delta;

		error = delta.rowwise().norm().maxCoeff();
		iter_num++;
		cout << "iter: " << iter_num << ", error: " << error << endl;
	}
}
// given Q_0,...,Q_m, fit by B-spline with control points P_0,...,P_n
void NURBSCurve::lspiafit(const MatrixXd & points, int n_cpts, int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->k = 4;
	const int m = points.rows() - 1;
	this->n = n_cpts - 1;
	const int dimension = points.cols();
	controlPw = MatrixXd::Zero(n_cpts, dimension);
	knots = VectorXd(n + k + 1);

	VectorXd params = parameterize(points);
	// initialize knot vector
	knots(0) = 0.0; knots(1) = 0.0; knots(2) = 0.0; knots(3) = 0.0;
	knots(n + 4) = 1.0; knots(n + 3) = 1.0; knots(n + 2) = 1.0; knots(n + 1) = 1.0;



	double d = 1.0 * (m + 1) / (n - 2);
	for (int j = 1; j <= n - 3; j++) {
		int i = j*d;
		double alpha = 1.0*j*d - i;
		knots(j + 3) = (1.0 - alpha)*params(i - 1) + alpha*params(i);
	}
	cout << "knots: " << knots.transpose() << endl;
	// initial control points P_0,...,P_n
	controlPw.row(0) = points.row(0);
	controlPw.row(n) = points.row(m);
	MatrixXd points_eval(points.rows(), points.cols());
	for (int j = 0; j < points.rows(); j++) {
		points_eval.row(j) = eval(params(j));
	}
	for (int iter = 0; iter < max_iter_num; iter++) {

		for (int i = 1; i < n; i++) {
			double sum1 = 0;
			VectorXd sum2 = VectorXd::Zero(dimension);

			for (int j = 0; j < params.size(); j++) {
				double blend = Basis(knots, params(j), i);
				sum1 += blend;
				VectorXd delta = points.row(j) - points_eval.row(j); // 拟合点处误差向量
				sum2 += blend * delta;
			}
			double factor = 0.0;
			if (abs(sum1) > 0.0001) {
				factor = 1.0 / sum1; // sum1为0时，对应点不更新
			}
			sum2 *= factor; // 控制点更新差向量
			controlPw.row(i) += sum2;

			double error = 0.0;
			// 控制点更新后，计算新的拟合位置和误差
			for (int j = 0; j < points.rows(); j++) {
				points_eval.row(j) = eval(params(j));
				error += (points.row(j) - points_eval.row(j)).norm();
			}
			error /= points.rows();
			//cout << "iter: " << iter + 1 << ", error: " << error << endl;
			if (error < eps) {
				break;
			}
		}

	}
	
}

void NURBSCurve::lspiafit(const MatrixXd & points, const VectorXd& params, int n_cpts, const VectorXd & knotvector, int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->k = 4;
	const int m = points.rows() - 1;
	this->n = n_cpts - 1;
	const int dimension = points.cols();
	controlPw = Eigen::MatrixXd::Zero(n_cpts, dimension);
	
	//VectorXd params = parameterize(points);

	knots = knotvector;

	controlPw.row(0) = points.row(0);
	controlPw.row(n) = points.row(m);
	MatrixXd points_eval(points.rows(), points.cols());
	for (int j = 0; j < points.rows(); j++) {
		points_eval.row(j) = eval(params(j));
	}
	for (int iter = 0; iter < max_iter_num; iter++) {
		
		for (int i = 1; i < n; i++) {
			double sum1 = 0;
			VectorXd sum2 = VectorXd::Zero(dimension);

			for (int j = 0; j < params.size(); j++) {
				double blend = Basis(knots, params(j), i);
				sum1 += blend;
				VectorXd delta = points.row(j) - points_eval.row(j); // 拟合点处误差向量
				sum2 += blend * delta;
			}
			double factor = 0.0;
			if (abs(sum1) > 0.0001) {
				factor = 1.0 / sum1; // sum1为0时，对应点不更新
			}
			sum2 *= factor; // 控制点更新差向量
			controlPw.row(i) += sum2;
			//cout << "controlpw " << i << " : " << controlPw.row(i) << endl;

			double error = 0.0;
			// 控制点更新后，计算新的拟合位置和误差
			for (int j = 0; j < points.rows(); j++) {
				points_eval.row(j) = eval(params(j));
				error += (points.row(j) - points_eval.row(j)).norm();
			}
			error /= points.rows();
			//cout << "iter: " << iter + 1 << ", error: " << error << endl;
			if (error < eps) {
				break;
			}
		}

	}
	
}
// 给定初始控制点
void NURBSCurve::lspiafit(const MatrixXd & points, const VectorXd & params, const MatrixXd & cpts, const VectorXd & knotvector, int max_iter_num, double eps)
{
	assert(points.rows() > 1 && points.cols() > 0);
	this->k = 4;
	const int m = points.rows() - 1;
	this->n = cpts.rows() - 1;
	const int dimension = points.cols();
	controlPw = cpts;
	

	knots = knotvector;

	controlPw.row(0) = points.row(0);
	controlPw.row(n) = points.row(m);
	MatrixXd points_eval(points.rows(), points.cols());
	for (int j = 0; j < points.rows(); j++) {
		points_eval.row(j) = eval(params(j));
	}
	for (int iter = 0; iter < max_iter_num; iter++) {

		for (int i = 1; i < n; i++) {
			double sum1 = 0;
			VectorXd sum2 = VectorXd::Zero(dimension);

			for (int j = 0; j < params.size(); j++) {
				double blend = Basis(knots, params(j), i);
				sum1 += blend;
				VectorXd delta = points.row(j) - points_eval.row(j); // 拟合点处误差向量
				sum2 += blend * delta;
			}
			double factor = 0.0;
			if (abs(sum1) > 0.0001) {
				factor = 1.0 / sum1; // sum1为0时，对应点不更新
			}
			sum2 *= factor; // 控制点更新差向量
			controlPw.row(i) += sum2;
			//cout << "controlpw " << i << " : " << controlPw.row(i) << endl;

			double error = 0.0;
			// 控制点更新后，计算新的拟合位置和误差
			for (int j = 0; j < points.rows(); j++) {
				points_eval.row(j) = eval(params(j));
				error += (points.row(j) - points_eval.row(j)).norm();
			}
			error /= points.rows();
			//cout << "iter: " << iter + 1 << ", error: " << error << endl;
			if (error < eps) {
				break;
			}
		}

	}
}

bool NURBSCurve::insert(double t)
{
	/*cout << "before inserting:------------------" << endl;
	cout << "controlPw:\n" << controlPw << endl;
	cout << "knots:\n" << knots << endl;*/
	assert(t >= knots(0) && t <= knots(n + k));
	int L = find_ind(t);
	// befor insert: p_0, ..., p_(L-k+1), p_(L-k+2),... ,		 p_L,...
	// after insert: p_0, ..., p_(L-k+1), p'_(L-k+2),..., p'_L,  p'_(L+1)=p_L,...   
	int start = (L-k+1>=0)?L-k+1:0;
	int end = (L<=n)?L:n;
	MatrixXd new_controlPw(controlPw.rows()+1,controlPw.cols());
	for(int i=0;i<new_controlPw.rows();i++){
		if(i<=start){
			new_controlPw.row(i) = controlPw.row(i);
		}else if(i<=end){
			double factor = (t-knots(i))/(knots(i+k-1)-knots(i));
			new_controlPw.row(i)=factor*controlPw.row(i)+(1.0-factor)*controlPw.row(i-1);
			
		}else{
			new_controlPw.row(i)=controlPw.row(i-1);

		}
	}

	VectorXd new_knots(knots.size()+1);
	for(int i=0;i<new_knots.size();i++){
		if(i<=L){
			new_knots(i)=knots(i);
		}else if(i==L+1){
			new_knots(i)=t;
		}else{
			new_knots(i)=knots(i-1);
		}
	}

	controlPw = new_controlPw;
	knots = new_knots;
	n+=1;
	/*cout << "after inserting:------------------" << endl;
	cout << "controlPw:\n" << controlPw << endl;
	cout << "knots:\n" << knots << endl;*/

	return true;
}


// draw controlpolygon
void NURBSCurve::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	
	// plot control points
	viewer.data().add_points(controlP, Eigen::RowVector3d(1, 1, 1));
	for (int i = 0; i <= n; i++) {
		viewer.data().add_label(controlP.row(i), std::to_string(i));
	}
	// plot control polygon
	for (int i = 0; i < n; i++){
		viewer.data().add_edges(
			controlP.row(i),
			controlP.row(i + 1),
			Eigen::RowVector3d(1, 1, 1));
	}
}

// draw NURBS surface
void NURBSCurve::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	double left = knots(k - 1);
	double right = knots(n + 1);
	const int num = (right - left) / resolution;
	double new_resolution = (right - left) / num;
	// plot NURBSCurve
	for (int i = 0; i < num; i++)
	{
		Eigen::MatrixXd P1 = eval(left + 1.0*i* new_resolution);
		Eigen::MatrixXd P2 = eval(left + 1.0*(i + 1)*new_resolution);

		if (isRational) {
			viewer.data().add_edges(
				P1.rowwise().hnormalized(),
				P2.rowwise().hnormalized(),
				Eigen::RowVector3d(1, 0, 0));
		}
		else {
			viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 0, 0));
		}

	}
}


// display by libigl
void NURBSCurve::draw(
	igl::opengl::glfw::Viewer& viewer, 
	bool showpolygon,bool showsurface,
	double resolution)
{
	if(isRational){
		controlP = controlPw.rowwise().hnormalized();
	}else{
		controlP = controlPw;
	}

	if(showpolygon){
		drawControlPolygon(viewer);
	}
	if(showsurface){
		drawSurface(viewer, resolution);
	}
	viewer.core.align_camera_center(controlP);
}





