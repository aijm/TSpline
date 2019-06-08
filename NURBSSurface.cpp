// This file is part of NURBS, a simple NURBS library.
// github repo: https://github.com/aijm/NURBS
// Copyright (C) 2018 Jiaming Ai <aichangeworld@gmail.com>

#include "NURBSSurface.h"


/*input format:
_order       : // _order(0):u direction; _order(1): v direction
_controlP: _controlP[i] represents u direction control point ,matrix (m+1) by 3 or 4
              v
              |
_controlP[n]: | P_0n P_1n ... P_mn
_controlP[i]: | ...
_controlP[1]: | P_01 P_11 ... P_m1
_controlP[0]: | P_00 P_10 ... P_m0
               ------------------------> u

_uknots   : u_0,u_1,...,u_(m+u_order)
_vknots   : v_0,v_1,...,v_(n+v_order)
_isRational:          */
NURBSSurface::NURBSSurface(
	VectorXi _order, 
	vector<MatrixXd> _controlP, 
	VectorXd _uknots, 
	VectorXd _vknots,
	bool _isRational)
{
	// pass and check the parameters of NURBS surface
	u_order = _order(0);
	v_order = _order(1);
	v_num = _controlP.size() - 1;
	assert(u_order >= 1 && v_order >= 1 && v_num >= 1);
	u_num = _controlP[0].rows() - 1;
	assert(_uknots.size() == u_num + u_order + 1 && _vknots.size() == v_num + v_order + 1);
	assert(u_num >= u_order - 1 && v_num >= v_order - 1);

	dimension = _controlP[0].cols(); // the dimension of control point 2 or 3 or 4(weighted)
	uknots = _uknots;
	vknots = _vknots;
	controlPw = _controlP;
	isRational = _isRational;

	if (isRational) { assert(dimension == 4); }
	else { assert(dimension == 3); }

}


// load
bool NURBSSurface::loadNURBS(string name){
	ifstream in(name);
	if(!in.is_open()){
		cout << "failed to open file: " + name << endl;
		return false;
	}
	char sep;
	dimension=3;
	in>> isRational;
	in>> u_num >> v_num;
	in>> u_order >> v_order;
	in>> dimension;
	controlPw = vector<MatrixXd>(v_num+1);
	for(int i=0;i<controlPw.size();i++){
		controlPw[i] = MatrixXd(u_num+1,dimension);
	}
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			for(int k=0;k<dimension;k++){
				in>> controlPw[i](j,k);
			}
			in>>sep;
		}
	}

	uknots = VectorXd(u_num+u_order+1);
	vknots = VectorXd(v_num+v_order+1);

	for(int i=0;i<uknots.size();i++){
		in>>uknots(i);
	}
	
	for(int i=0;i<vknots.size();i++){
		in>>vknots(i);
	}
	return true;
}

// save 
bool NURBSSurface::saveNURBS(string name){
	if(controlPw.size()==0){
		cout<< "nothing to save!"<<endl;
		return false;
	}
	if(isRational){
		name+=".cptw";
	}else{
		name+=".cpt";
	}
	ofstream out(name);
	if(!out.is_open()){
		cout << "failed to open or create file: " + name << endl;
		return false;
	}
	IOFormat outputFmt(4, 0, " ", " ", "", ",");	
	out<< isRational<<endl;
	out<< u_num << " " <<v_num<<endl;
	out<< u_order << " "<< v_order<<endl;
	out<< controlPw[0].cols()<<endl;
	
	for(int i=0;i<controlPw.size();i++){
		for(int j=0;j<=u_num;j++){
			out<<controlPw[i].row(j).format(outputFmt);
		}
		out<<endl;
	}
	out<<uknots.transpose()<<endl;
	out<<vknots.transpose();
	return true;
}

bool NURBSSurface::saveAsObj(string filename, double resolution)
{
	filename += ".obj";
	ofstream out(filename);
	if (!out.is_open()) {
		cout << "can't open file: " << filename << endl;
		return false;
	}

	double u_low = uknots(u_order - 1);
	double u_high = uknots(u_num + 1);
	const int uspan = (u_high - u_low) / resolution;
	double u_resolution = (u_high - u_low) / uspan;

	double v_low = vknots(v_order - 1);
	double v_high = vknots(v_num + 1);
	const int vspan = (v_high - v_low) / resolution;
	double v_resolution = (v_high - v_low) / vspan;

	mesh_V = MatrixXd((uspan + 1)*(vspan + 1), 3);
	mesh_F = MatrixXi(2 * uspan*vspan, 3);
	// discretize NURBS Surface into triangular mesh(V,F) in libigl mesh structure
	// calculate mesh_V
	for (int j = 0; j <= vspan; j++)
		for (int i = 0; i <= uspan; i++)
		{
			RowVectorXd curvePoint = eval(u_low + i*u_resolution, v_low + j*v_resolution).row(0);
			//cout << "curvepoint: " << curvePoint << endl;
			if (isRational) { mesh_V.row(j*(uspan + 1) + i) = curvePoint.hnormalized(); }
			else { mesh_V.row(j*(uspan + 1) + i) = curvePoint; }
		}

	for (int j = 0; j<vspan; j++)
		for (int i = 0; i < uspan; i++)
		{
			int V_index = j*(uspan + 1) + i;
			int F_index = 2 * j*uspan + 2 * i;
			mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
			mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
		}

	igl::writeOBJ(filename, mesh_V, mesh_F);
	return true;
}
	

// find the knot interval of t by binary searching
int NURBSSurface::find_ind(double t, int k, int n, const VectorXd& knots)
{
	if (t == knots(n + 1)) return n;
	int low = 0;
	int high = n + k;
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
// calculate coordinate of curve point with parameter u & v
MatrixXd NURBSSurface::eval(double u, double v) const
{
	MatrixXd res = MatrixXd::Zero(1, dimension);
	// 只需计算 4*4个点
	int uid = FindSpan(uknots, u);
	int vid = FindSpan(vknots, v);

	for (int i = uid - 1; i <= uid + 2; i++) {
		for (int j = vid - 1; j <= vid + 2; j++) {
			double blend = NURBSCurve::Basis(uknots, u, i - 2)*NURBSCurve::Basis(vknots, v, j - 2);
			res += controlPw[j - 2].row(i - 2) * blend;
		}
	}
	return res;
	////Calculating the Control Points of U-direction Isoparametric Line
	//MatrixXd v_controlP(v_num + 1, dimension); 
	//for (int i = 0; i <= v_num; i++)
	//{
	//	v_controlP.row(i) = eval(u, controlPw[i], uknots);
	//}
	//return eval(v, v_controlP, vknots); // Calculating the coordinate of parameter v
}

MatrixXd NURBSSurface::eval(
	double t, 
	const MatrixXd &_controlP, 
	const VectorXd &knots)
{
	
	int n = _controlP.rows() - 1;
	int k = knots.size() - _controlP.rows();
	/*if (t<knots(k - 1) || t>knots(n + 1)) {*/
		//cout << "t: " << t << ", left: " << knots(k - 1) << ", right: " << knots(n + 1) << endl;
	/*}*/
	assert(t>=knots(k-1) && t<=knots(n+1));
	// find the knot interval of t by binary searching
	int L = find_ind(t, k, n, knots); //[t_L,t_(L+1)] 

	// P_(L-k+1),..,P_L control the interval [t_L,t_(L+1)]
	MatrixXd temp = _controlP.block(L - k + 1, 0, k, _controlP.cols());

	for (int r = 1; r <= k - 1; r++)
		for (int i = L - k + 1 + r; i <= L; i++) 
		{
			double factor = (t - knots(i)) / (knots(i + k - r) - knots(i));
			int start = i - (L - k + 1 + r);
			temp.row(start) = (1.0 - factor)*temp.row(start) + factor*temp.row(start + 1);
		}
	MatrixXd curvePoint = temp.row(0);
	return curvePoint;
}

// knot insertion
bool NURBSSurface::insert(double s, char dir){
	assert(dir=='u' || dir=='v');
	if(dir=='u'){
		for(int i=0;i<controlPw.size();i++){
			NURBSCurve nurbs(u_num,u_order,controlPw[i],uknots,isRational);
			nurbs.insert(s);
			controlPw[i]=nurbs.controlPw;
			if (i == controlPw.size() - 1) {
				uknots = nurbs.knots;
				u_num = nurbs.n;
			}
		}
		return true;
	}else if(dir=='v'){
		vector<MatrixXd> new_controlPw(controlPw.size()+1);
		for (int i = 0; i < new_controlPw.size(); i++) {
			new_controlPw[i] = MatrixXd(u_num + 1, controlPw[0].cols());
		}
		MatrixXd v_controlPw(v_num+1,controlPw[0].cols());
		for(int i=0;i<=u_num;i++){
			for(int j=0;j<=v_num;j++){
				v_controlPw.row(j) = controlPw[j].row(i);
			}
			NURBSCurve nurbs(v_num,v_order,v_controlPw,vknots,isRational);
			nurbs.insert(s);
			
			for(int k=0;k<new_controlPw.size();k++){
				new_controlPw[k].row(i) = nurbs.controlPw.row(k);
			}
			if (i == u_num) {
				vknots = nurbs.knots;
				v_num = nurbs.n;
			}
		}
		controlPw = new_controlPw;
		return true;

	}else{
		cout<< "please input dir as u or v!"<<endl;
		return false;
	}
}

// kont insertion
bool NURBSSurface::insert(double s, double t){
	if(insert(s,'u') && insert(t,'v')){
		return true;
	}
	return false;
}

// draw controlpolygon
void NURBSSurface::drawControlPolygon(igl::opengl::glfw::Viewer &viewer){
	
	//plot control points and control polygon
	
	//plot control points
	for (int i = 0; i < controlP.size(); i++)
	{
		viewer.data().add_points(
			controlP[i],
			Eigen::RowVector3d(1, 1, 1));
	}
	// plot control polygon
	for (int j = 0; j <= v_num; j++)
		for (int i = 0; i <= u_num; i++)
		{
			if (i != u_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j].row(i + 1),
					Eigen::RowVector3d(1, 1, 1));
			}
			if (j != v_num) {
				viewer.data().add_edges(
					controlP[j].row(i),
					controlP[j + 1].row(i),
					Eigen::RowVector3d(1, 1, 1));
			}
		}
}

// draw NURBS surface
void NURBSSurface::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution){
	// cut apart the parameter domain
	double u_low = uknots(u_order - 1);
	double u_high = uknots(u_num + 1);
	const int uspan = (u_high - u_low) / resolution;
	double u_resolution = (u_high - u_low) / uspan;

	double v_low = vknots(v_order - 1);
	double v_high = vknots(v_num + 1);
	const int vspan = (v_high - v_low) / resolution;
	double v_resolution = (v_high - v_low) / vspan;

	mesh_V = MatrixXd((uspan + 1)*(vspan + 1), 3);
	mesh_F = MatrixXi(2 * uspan*vspan, 3);

	VectorXd HH = VectorXd::Zero(mesh_V.rows());
	VectorXd KK = VectorXd::Zero(mesh_V.rows());
	VectorXd PV11 = VectorXd::Zero(mesh_V.rows());
	VectorXd PV22 = VectorXd::Zero(mesh_V.rows());
	// discretize NURBS Surface into triangular mesh(V,F) in libigl mesh structure
	// calculate mesh_V
	for (int j = 0; j <= vspan; j++)
		for (int i = 0; i <= uspan; i++)
		{
			double u = u_low + i*u_resolution;
			double v = v_low + j*v_resolution;
			RowVectorXd curvePoint = eval(u, v).row(0);

			RowVector3d du, dv, d2u, d2v, duv;
			derivative(u, v, du, dv, d2u, d2v, duv);
			/*cout << "d2u: " << d2u << endl;
			double h = 0.001;*/
			//RowVector3d d2u_1 = (eval(u + h, v) - 2 * eval(u, v) + eval(u - h, v)) / (h*h);
			
			//cout << "d2u_1: " << d2u_1 << endl;
			RowVector3d normal = du.cross(dv);
			normal.normalize();
			/*viewer.data().add_edges(curvePoint, curvePoint + 0.005*du.normalized(), Eigen::RowVector3d(1, 0, 0));
			viewer.data().add_edges(curvePoint, curvePoint + 0.005*dv.normalized(), Eigen::RowVector3d(0, 1, 0));
			viewer.data().add_edges(curvePoint, curvePoint + 0.005*normal, Eigen::RowVector3d(0, 0, 1));
			*/

			int id = j*(uspan + 1) + i;
			KK(id) = guassian_curvature(u, v);

			if (isRational) { mesh_V.row(id) = curvePoint.hnormalized(); }
			else { mesh_V.row(id) = curvePoint; }
		}

	for (int j = 0; j<vspan; j++)
		for (int i = 0; i < uspan; i++)
		{
			int V_index = j*(uspan + 1) + i;
			int F_index = 2 * j*uspan + 2 * i;
			mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
			mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
		}

	VectorXd K;
	// Compute integral of Gaussian curvature
	igl::gaussian_curvature(mesh_V, mesh_F, K);
	// Compute mass matrix
	SparseMatrix<double> M, Minv;
	igl::massmatrix(mesh_V, mesh_F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	igl::invert_diag(M, Minv);
	// Divide by area to get integral average
	K = (Minv*K).eval();
	KK = KK.cwiseAbs();
	for (int i = 0; i < KK.size(); i++) {
		KK(i) = log(KK(i) + 1);
	}

	Eigen::VectorXd H;

							 // compute curvatrue directions via quadric fitting
	Eigen::MatrixXd PD1, PD2;
	Eigen::VectorXd PV1, PV2;
	igl::principal_curvature(mesh_V, mesh_F, PD1, PD2, PV1, PV2);
	// mean curvature
	H = 0.5*(PV1 + PV2);

	for (int i = 0; i < K.rows(); i++) {
		cout << K(i) << ", " << KK(i) << endl;
	}
	

	viewer.data().set_mesh(mesh_V, mesh_F);

	// Compute pseudocolor
	Eigen::MatrixXd C;
	igl::parula(KK, true, C);
	viewer.data().set_colors(C);
}


void NURBSSurface::draw(
	igl::opengl::glfw::Viewer &viewer, 
	bool showpolygon,bool showsurface,
	double resolution){

	if(controlP.size()!=controlPw.size()){
		controlP = vector<MatrixXd>(controlPw.size());
	}
	if(isRational){
		for(int i=0;i<controlP.size();i++){
			controlP[i] = controlPw[i].rowwise().hnormalized();
		}
	}else{
		controlP = controlPw;
	}
	
	if(showpolygon){
		drawControlPolygon(viewer);
		viewer.core.align_camera_center(controlP[controlP.size()/2]);
	}
	if(showsurface){
		drawSurface(viewer,resolution);
		viewer.core.align_camera_center(mesh_V,mesh_F);
	}

}

// surface skinning: order is the same for every curves
void NURBSSurface::skinning(const vector<NURBSCurve> &curves,igl::opengl::glfw::Viewer &viewer){
	assert(curves.size()>1);
	isRational = curves[0].isRational;
	// count knots of every curve
	vector<map<double,int>> curve_knots(curves.size());
	for(int i=0;i<curves.size();i++){
		for(int j=0;j<curves[i].knots.size();j++){
			if(curve_knots[i].find(curves[i].knots(j))==curve_knots[i].end()){
				curve_knots[i][curves[i].knots(j)] = 1;
			}else{
				curve_knots[i][curves[i].knots(j)] += 1;
			}
		}
	}
	// compute merge knots 
	map<double, int> merge_knots;
	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = curve_knots[i].begin(); it != curve_knots[i].end(); it++) {
			if (merge_knots.find(it->first) == merge_knots.end()) {
				merge_knots[it->first] = it->second;
			}
			else if (merge_knots[it->first] < it->second) {
				merge_knots[it->first] = it->second;
			}
		}
	}
	//cout << "merge_knots:\n" << endl;
	/*for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
		cout << it->first << ", " << it->second << endl;
	}*/

	// change every curve to 3D and rational
	vector<NURBSCurve> new_curves = curves;

	// for every curve, insert knots to make curve knots compatible
	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
			int insert_num = 0;
			if (curve_knots[i].find(it->first) == curve_knots[i].end()) {
				insert_num = it->second;
			}
			else {
				insert_num = it->second - curve_knots[i][it->first];
			}
			// insert
			for (int j = 0; j < insert_num; j++) {
				new_curves[i].insert(it->first);
			}
		}
	}
	
	for (int i = 0; i < new_curves.size(); i++) {
		//cout << "curve " << i << ": " << new_curves[i].knots.transpose() << endl;
		/*if (i == 3 || i == 4) {
			viewer.data().add_points(new_curves[i].controlPw.rowwise().hnormalized(), RowVector3d(0, 1, 0));
		}*/
	}
	

	vknots = new_curves[0].knots;
	//cout << "after insert knots: " << vknots.transpose() << endl;
	// compute u-direction knot vector
	VectorXd curves_param = VectorXd::Zero(new_curves.size());
	v_num = new_curves[0].n;
	v_order = new_curves[0].k;
	dimension = new_curves[0].controlPw.cols();
	//u_num = new_curves.size() - 1;
	u_order = 4;
	int curve_num = curves.size() - 1;
	u_num = curve_num + 2;

	for (int i = 0; i <= v_num; i++) {
		MatrixXd u_cpts(curve_num + 1, dimension);
		for (int j = 0; j <= curve_num; j++) {
			u_cpts.row(j) = new_curves[j].controlPw.row(i);
		}
		curves_param += NURBSCurve::parameterize(u_cpts);
	}
	curves_param /= (v_num + 1);
	uknots = VectorXd::Zero(curve_num + 7); // u_0,u_0,u_0, u_0,...,u_K, u_K,u_K,u_K

	uknots.block(3, 0, curve_num + 1, 1) = curves_param;
	uknots(curve_num + 6) = 1.0; uknots(curve_num + 5) = 1.0; uknots(curve_num + 4) = 1.0;

	//uknots << 0, 0, 0, 0, 2.0/12, 4.0/12, 5.0/12, 7.0/12, 8.0/12, 11.0/12, 1, 1, 1, 1;

	//cout << "uknots: " << uknots.transpose() << endl;
	// interpolate points to get the NURBS skinning surface control points
	controlPw = vector<MatrixXd>(v_num + 1);
	for (int i = 0; i < controlPw.size(); i++) {
		controlPw[i] = MatrixXd(u_num + 1, dimension);
	}

	for (int i = 0; i <= v_num; i++) {
		MatrixXd u_ctps = MatrixXd::Zero(curve_num + 1, dimension);
		for (int j = 0; j <= curve_num; j++) {
			u_ctps.row(j) = new_curves[j].controlPw.row(i);
		}
		NURBSCurve nurbs;
		
		nurbs.interpolate(u_ctps, uknots);

		//nurbs.isRational = false;
		/*if (i==3) {
			viewer.data().add_points(u_ctps.rowwise().hnormalized(), RowVector3d(0, 1, 0));
			cout << "points:\n" << u_ctps << endl;
			cout << "knots: " << nurbs.knots.transpose() << endl;
			nurbs.draw(viewer, true, true);
			
		}*/
		
		//cout << "dim: " << nurbs.controlPw.cols() << endl;
		controlPw[i] = nurbs.controlPw;
		/*if (i == 3) {
			viewer.data().add_points(controlPw[i].rowwise().hnormalized(), RowVector3d(1, 0, 0));
		}
		else if (i == 4) {
			viewer.data().add_points(controlPw[i].rowwise().hnormalized(), RowVector3d(0, 1, 0));
		}*/
	}
	//dimension = 3;
	cout << "number of control points: " << controlPw.size()*controlPw[0].rows() << endl;
	cout << "skinning finished!" << endl;
}

void NURBSSurface::skinning(const vector<NURBSCurve>& curves, const VectorXd & curves_param, igl::opengl::glfw::Viewer & viewer)
{
	assert(curves.size()>1);
	isRational = curves[0].isRational;
	// count knots of every curve
	vector<map<double, int>> curve_knots(curves.size());
	for (int i = 0; i<curves.size(); i++) {
		for (int j = 0; j<curves[i].knots.size(); j++) {
			if (curve_knots[i].find(curves[i].knots(j)) == curve_knots[i].end()) {
				curve_knots[i][curves[i].knots(j)] = 1;
			}
			else {
				curve_knots[i][curves[i].knots(j)] += 1;
			}
		}
	}
	// compute merge knots 
	map<double, int> merge_knots;
	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = curve_knots[i].begin(); it != curve_knots[i].end(); it++) {
			if (merge_knots.find(it->first) == merge_knots.end()) {
				merge_knots[it->first] = it->second;
			}
			else if (merge_knots[it->first] < it->second) {
				merge_knots[it->first] = it->second;
			}
		}
	}
	//cout << "merge_knots:\n" << endl;
	/*for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
	cout << it->first << ", " << it->second << endl;
	}*/

	// change every curve to 3D and rational
	vector<NURBSCurve> new_curves = curves;

	// for every curve, insert knots to make curve knots compatible
	for (int i = 0; i < curve_knots.size(); i++) {
		for (auto it = merge_knots.begin(); it != merge_knots.end(); it++) {
			int insert_num = 0;
			if (curve_knots[i].find(it->first) == curve_knots[i].end()) {
				insert_num = it->second;
			}
			else {
				insert_num = it->second - curve_knots[i][it->first];
			}
			// insert
			for (int j = 0; j < insert_num; j++) {
				new_curves[i].insert(it->first);
			}
		}
	}

	for (int i = 0; i < new_curves.size(); i++) {
		//cout << "curve " << i << ": " << new_curves[i].knots.transpose() << endl;
		/*if (i == 3 || i == 4) {
		viewer.data().add_points(new_curves[i].controlPw.rowwise().hnormalized(), RowVector3d(0, 1, 0));
		}*/
	}


	vknots = new_curves[0].knots;
	//cout << "after insert knots: " << vknots.transpose() << endl;
	// compute u-direction knot vector
	//VectorXd curves_param = VectorXd::Zero(new_curves.size());
	v_num = new_curves[0].n;
	v_order = new_curves[0].k;
	dimension = new_curves[0].controlPw.cols();
	//u_num = new_curves.size() - 1;
	u_order = 4;
	int curve_num = curves.size() - 1;
	u_num = curve_num + 2;

	uknots = VectorXd::Zero(curve_num + 7); // u_0,u_0,u_0, u_0,...,u_K, u_K,u_K,u_K

	uknots.block(3, 0, curve_num + 1, 1) = curves_param;
	uknots(curve_num + 6) = 1.0; uknots(curve_num + 5) = 1.0; uknots(curve_num + 4) = 1.0;

	//uknots << 0, 0, 0, 0, 2.0/12, 4.0/12, 5.0/12, 7.0/12, 8.0/12, 11.0/12, 1, 1, 1, 1;

	//cout << "uknots: " << uknots.transpose() << endl;
	// interpolate points to get the NURBS skinning surface control points
	controlPw = vector<MatrixXd>(v_num + 1);
	for (int i = 0; i < controlPw.size(); i++) {
		controlPw[i] = MatrixXd(u_num + 1, dimension);
	}

	for (int i = 0; i <= v_num; i++) {
		MatrixXd u_ctps = MatrixXd::Zero(curve_num + 1, dimension);
		for (int j = 0; j <= curve_num; j++) {
			u_ctps.row(j) = new_curves[j].controlPw.row(i);
		}
		NURBSCurve nurbs;

		nurbs.interpolate(u_ctps, uknots);

		//nurbs.isRational = false;
		/*if (i==3) {
		viewer.data().add_points(u_ctps.rowwise().hnormalized(), RowVector3d(0, 1, 0));
		cout << "points:\n" << u_ctps << endl;
		cout << "knots: " << nurbs.knots.transpose() << endl;
		nurbs.draw(viewer, true, true);

		}*/

		//cout << "dim: " << nurbs.controlPw.cols() << endl;
		controlPw[i] = nurbs.controlPw;
		/*if (i == 3) {
		viewer.data().add_points(controlPw[i].rowwise().hnormalized(), RowVector3d(1, 0, 0));
		}
		else if (i == 4) {
		viewer.data().add_points(controlPw[i].rowwise().hnormalized(), RowVector3d(0, 1, 0));
		}*/
	}
	//dimension = 3;
	cout << "number of control points: " << controlPw.size()*controlPw[0].rows() << endl;
	cout << "skinning finished!" << endl;
}

double NURBSSurface::mean_curvature(double u, double v) const
{
	double k1, k2;
	curvature(u, v, k1, k2);
	return (k1 + k2) / 2;
}

double NURBSSurface::guassian_curvature(double u, double v) const
{
	double k1, k2;
	curvature(u, v, k1, k2);
	return k1*k2;
}

void NURBSSurface::curvature(double u, double v, double & k1, double & k2) const
{
	RowVector3d du = RowVector3d::Zero();
	RowVector3d dv = RowVector3d::Zero();
	RowVector3d d2u = RowVector3d::Zero();
	RowVector3d d2v = RowVector3d::Zero();
	RowVector3d duv = RowVector3d::Zero();
	derivative(u, v, du, dv, d2u, d2v, duv);
	RowVector3d normal = du.cross(dv);
	normal.normalize();
	double E = du.dot(du);
	double F = du.dot(dv);
	double G = dv.dot(dv);
	double L = d2u.dot(normal);
	double M = duv.dot(normal);
	double N = d2v.dot(normal);
	double b = G*L - 2 * F*M + E*N;
	double a = E*G - F*F;
	double c = L*N - M*M;
	double temp = sqrt(b*b - 4 * a*c);
	k1 = 0.5*(b - temp) / a;
	k2 = 0.5*(b + temp) / a;
}

void NURBSSurface::derivative(double u, double v, RowVector3d & du, RowVector3d & dv, RowVector3d & d2u, RowVector3d & d2v, RowVector3d & duv) const
{
	du = RowVector3d::Zero();
	dv = RowVector3d::Zero();
	d2u = RowVector3d::Zero();
	d2v = RowVector3d::Zero();
	duv = RowVector3d::Zero();
	
	double h = 0.001;
	if (u < h) u += h;
	if (u > 1 - h) u -= h;
	if (v < h) v += h;
	if (v > 1 - h) v -= h;

	RowVector3d P00 = eval(u - h, v - h);
	RowVector3d P10 = eval(u, v - h);
	RowVector3d P20 = eval(u + h, v - h);
	RowVector3d P01 = eval(u - h, v);
	RowVector3d P11 = eval(u, v);
	RowVector3d P21 = eval(u + h, v);

	RowVector3d P02 = eval(u - h, v + h);
	RowVector3d P12 = eval(u, v + h);
	RowVector3d P22 = eval(u + h, v + h);

	du = (P21 - P01) / 2.0 / h;
	dv = (P12 - P10) / 2.0 / h;
	d2u = (P01 - 2 * P11 + P21) / h / h;
	d2v = (P10 - 2 * P11 + P12) / h / h;
	duv = (P22 - P02 - P20 + P00) / 4 / h / h;

	/*else {
		int uid = FindSpan(uknots, u);
		int vid = FindSpan(uknots, v);

		for (int i = uid - 1; i <= uid + 2; i++) {
			for (int j = uid - 1; j <= uid + 2; j++) {
				VectorXd der_u = NURBSCurve::DersBasis(uknots, u, i - 2);
				VectorXd der_v = NURBSCurve::DersBasis(vknots, v, j - 2);

				du += controlPw[j - 2].row(i - 2) * der_u(1)*der_v(0);
				dv += controlPw[j - 2].row(i - 2) * der_u(0)*der_v(1);
				d2u += controlPw[j - 2].row(i - 2) * der_u(2)*der_v(0);
				d2v += controlPw[j - 2].row(i - 2) * der_u(0)*der_v(2);
				duv += controlPw[j - 2].row(i - 2) * der_u(1)*der_v(1);

			}
		}
	}*/
	
	
}

int NURBSSurface::FindSpan(const Eigen::MatrixXd & knots, double t, int p)
{
	const int n = knots.size() - p - 2;
	if (t == knots(n + 1)) return n;
	int low = p;
	int high = n + 1;
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
