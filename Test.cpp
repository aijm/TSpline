#include "Test.h"

clock_t Test::begin;
clock_t Test::end;

void Test::test_Nurbs_curvature()
{
	NURBSSurface surface;
	surface.loadNURBS("../out/nurbs/venus_front.cpt");
	surface.draw(Window::viewer, false, true, 0.01);
	Window w;
	w.launch();
}

void Test::test_TsplineSimplify()
{
	vector<int> a{ 8,9,7,4,3,5 };
	vector<int> index(a.size());
	iota(index.begin(), index.end(), 0);
	sort(index.begin(), index.end(), [a](int id1, int id2) {return a[id1] > a[id2]; });
	for (int i = 0; i < index.size(); i++) {
		cout << "id: " << index[i] << "; value: " << a[index[i]] << endl;
	}
	NURBSSurface surface;
	surface.loadNURBS("../out/nurbs/venus_front.cpt");
	//surface.draw(Window::viewer, false, true);
	/*Window w;
	w.launch();*/
	Mesh3d tspline;   
	t_mesh::TsplineSimplify(surface, tspline,20, 5e-3);
	tspline.saveMesh("../out/tspline/fuck");
	MeshRender render(&tspline, false, true, true);
	render.launch();
}

void Test::test_generate_curves()
{
	BsplineVolume volume;
	volume.readVolume("../out/volume/venus_bspline.txt");
	const int curves_num = volume.knot_vector[2].size() - 6;
	vector<NURBSCurve> curves(curves_num);

	uniform_int_distribution<unsigned> u(15, 30);
	default_random_engine e;
	default_random_engine e1;
	uniform_real_distribution<double> u1(0, 1);

	for (int i = 0; i < curves_num; i++) {
		int sampleNum = u(e);
		cout << i << " sample " << sampleNum << " points" << endl;
		Eigen::MatrixXd points(sampleNum+2, 3);
		double z = volume.knot_vector[2](i + 3);

		vector<double> x_coors;
		x_coors.push_back(0);
		x_coors.push_back(1);
		
		for (int j = 0; j < sampleNum; j++) {
			x_coors.push_back(u1(e1));
		}
		sort(x_coors.begin(), x_coors.end());
		for (int j = 0; j < x_coors.size(); j++) {
			points.row(j) = volume.eval(x_coors[j], 0, z).toVectorXd();
		}
		curves[i].interpolate(points);
		curves[i].saveNURBS("../out/nurbs/venus_curve_" + to_string(i));
		curves[i].draw(Window::viewer, false, true, 0.001);
	}


	//Skinning* method = new PiaMethod(curves, 1000);
	//Skinning* method = new NasriMethod(curves);
	//Skinning* method = new OptMethod(curves);
	//Skinning* method = new PiaMinJaeMethod(curves, 1000);
	//Skinning* method = new PiaNasriMethod(curves, 1000);
	Skinning* method = new MinJaeMethod(curves, 40, 20);

	method->setViewer(&Window::viewer);
	method->calculate();
	Mesh3d* mesh = &(method->tspline);
	cout << "num of nodes: " << mesh->get_num() << endl;

	mesh->saveMesh("../out/tspline/venus_skinning");
	MeshRender render(mesh, false, false, true, 0.001);
	render.launch();

}

void Test::test_generate_curves1()
{
	/*BsplineVolume volume;
	volume.readVolume("../out/volume/venus_bspline.txt");*/

	NURBSSurface surface;

	
	surface.loadNURBS("../out/nurbs/venus_front.cpt");
	surface.saveAsObj("../out/OBJ/venus_surface", 0.01);
	//surface.draw(Window::viewer, false, true, 0.01);
	////surface.saveNURBS("../out/nurbs/venus_front");
	//Window w;
	//w.launch();

	const double resolution = 0.1;
	const int curves_num = 1.0 / resolution + 1;
	vector<NURBSCurve> curves(curves_num);

	uniform_int_distribution<unsigned> u(15, 30);
	default_random_engine e;
	default_random_engine e1;
	uniform_real_distribution<double> u1(0, 1);

	for (int i = 0; i < curves_num; i++) {
		int sampleNum = u(e);
		cout << i << " sample " << sampleNum << " points" << endl;
		Eigen::MatrixXd points(sampleNum + 2, 3);
		double z = i*resolution;

		vector<double> x_coors;
		x_coors.push_back(0);
		x_coors.push_back(1);

		for (int j = 0; j < sampleNum; j++) {
			double x = 1.0*(j+1) / (sampleNum+1);
			//double x = u1(e1);
			//cout << "x: " << x << endl;
			x_coors.push_back(x);
		}
		sort(x_coors.begin(), x_coors.end());
		for (int j = 0; j < x_coors.size(); j++) {
			points.row(j) = surface.eval(x_coors[j], z);
		}
		//Window::viewer.data().add_points(points, blue);
		curves[i].interpolate(points);
		curves[i].saveNURBS("../out/nurbs/venus_curve_1_" + to_string(i));
		curves[i].draw(Window::viewer, false, true,0.001);
	}

	/*NURBSSurface skin;
	skin.skinning(curves, Window::viewer);
	skin.draw(Window::viewer, false, true);
	skin.saveAsObj("../out/OBJ/venus_nurbs_skin");
	Window w;
	w.launch();*/
	//Skinning* method = new PiaMethod(curves, 100);
	//Skinning* method = new NasriMethod(curves);
	//Skinning* method = new OptMethod(curves);
	PiaMinJaeMethod* method = new PiaMinJaeMethod(curves, 100);
	//PiaNasriMethod* method = new PiaNasriMethod(curves, 100);
	//Skinning* method = new MinJaeMethod(curves, 40, 10);

	method->setViewer(&Window::viewer);
	MatrixXd helper_points;
	loadpoints("../out/points/helper_points.dat", helper_points);
	//method->set_helper_points(helper_points);

	method->calculate();
	Mesh3d* mesh = &(method->tspline);
	cout << "num of nodes: " << mesh->get_num() << endl;

	mesh->saveMesh("../out/tspline/venus_skinning1");
	//mesh->saveAsObj("../out/OBJ/venus_PiaMinJae_sample_400", 0.01);
	MeshRender render(mesh,false,true,true,0.01);
	render.launch();
}

void Test::test_generate_curves2()
{
	NURBSSurface surface;


	surface.loadNURBS("../out/nurbs/Bsurface_standard.cpt");
	//surface.draw(Window::viewer, false, true);
	const double resolution = 0.05;
	const int curves_num = 1.0 / resolution;
	cout << curves_num << endl;
	vector<NURBSCurve> curves(curves_num + 1);

	for (int i = 0; i <= curves_num; i++) {
		double t = 1.0*i / curves_num;
		surface.get_isoparam_curve(curves[i], t, 'u');
		curves[i].saveNURBS("../out/nurbs/Bsurface_curve_" + to_string(i));
		curves[i].draw(Window::viewer, false, true, 0.001);
	}

	/*Window w;
	w.launch();*/
	//Skinning* method = new PiaMethod(curves, 100);
	//Skinning* method = new NasriMethod(curves);
	//Skinning* method = new OptMethod(curves);
	PiaMinJaeMethod* method = new PiaMinJaeMethod(curves, 100);
	//PiaNasriMethod* method = new PiaNasriMethod(curves, 100);
	//Skinning* method = new MinJaeMethod(curves, 40, 10);

	method->setViewer(&Window::viewer);

	method->calculate();
	Mesh3d* mesh = &(method->tspline);
	cout << "num of nodes: " << mesh->get_num() << endl;

	mesh->saveMesh("../out/tspline/Bsurface_skinning1");
	mesh->saveAsObj("../out/OBJ/Bsurface_PiaMinJae_20", 0.01);
	MeshRender render(mesh, false, false, true, 0.01);
	render.launch();
}

void Test::load_nurbs_surface(NURBSSurface & surface, string filename)
{
	// 加载另外一种文件格式的nurbs surface
	ifstream in(filename);
	if (!in.is_open()) {
		cout << "error: can't open file: " << filename << endl;
		return;
	}
	
	surface.isRational = false;
	surface.u_order = 4;
	surface.v_order = 4;
	surface.dimension = 3;

	string line;
	getline(in, line);
	getline(in, line);
	stringstream ss0(line);
	ss0 >> surface.v_num >> surface.u_num;
	--surface.v_num;
	--surface.u_num;

	getline(in, line);
	surface.controlPw.resize(surface.v_num + 1);
	for (int i = 0; i <= surface.v_num; i++) {
		surface.controlPw[i] = MatrixXd::Zero(surface.u_num + 1, 3);
		for (int j = 0; j <= surface.u_num; j++) {
			getline(in, line);
			istringstream ss(line);

			ss >> surface.controlPw[i](j, 0) >> surface.controlPw[i](j, 1) >> surface.controlPw[i](j, 2);
		
			//cout << surface.controlPw[i](j, 0) << ", " << surface.controlPw[i](j, 1) << ", " << surface.controlPw[i](j, 2) << endl;
		}
	}
	//cout << "***********************************************" << endl;
	getline(in, line);
	getline(in, line);
	istringstream ss1(line);
	surface.vknots = VectorXd::Zero(surface.v_num + 5);
	for (int i = 0; i < surface.vknots.size(); i++) {
		ss1 >> surface.vknots(i);
	}
	//cout << "vkonts: \n" << surface.vknots << endl;

	getline(in, line);
	getline(in, line);
	istringstream ss2(line);
	surface.uknots = VectorXd::Zero(surface.u_num + 5);
	for (int i = 0; i < surface.uknots.size(); i++) {
		ss2 >> surface.uknots(i);
	}
	//cout << "ukonts: \n" << surface.uknots << endl;
	in.close();
}

void Test::test_load_nurbs_surface()
{
	string filename = "../out/nurbs/Bsurface.cpt";
	NURBSSurface surface;
	load_nurbs_surface(surface, filename);
	surface.draw(Window::viewer, false, true);
	surface.saveNURBS("../out/nurbs/Bsurface_standard");
	Window w;
	w.launch();
}

/**
   棋子蒙皮
*/
void Test::test_chess()
{
	
	vector<MatrixXd> curve_points(13);
	// 读取文件，得到散点
	for (int i = 0; i < curve_points.size(); i++) {
		string filename = "../out/curves/selected_point" + to_string(i + 1) + ".txt";
		ifstream in(filename);
		if (!in.is_open()) {
			cout << "error : can't open file: " << filename << endl;
			return;
		}
		string line;
		char c; // 读取'v'
		while (getline(in, line)) {
			curve_points[i].conservativeResize(curve_points[i].rows() + 1, 3);
			int lastid = curve_points[i].rows() - 1;
			stringstream ss(line);
			ss >> c >> curve_points[i](lastid, 0) >> curve_points[i](lastid, 1) >> curve_points[i](lastid, 2);
		}
		in.close();
	}
	
	// 插值得到B样条曲线
	vector<NURBSCurve> nurbs_curves(curve_points.size());
	for (int i = 0; i < curve_points.size(); i++) {
		nurbs_curves[i].interpolate(curve_points[i]);
		nurbs_curves[i].draw(Window::viewer, false, true);
		Window::viewer.data().add_points(curve_points[i], green);
	}

	Skinning* method = new PiaMinJaeMethod(nurbs_curves, 20);
	method->setViewer(&Window::viewer);
	method->calculate();
	MeshRender render(&method->tspline, false, true, true);
	render.launch();

}

/**
 从B样条体中生成一组B样条曲面
 用于后续简化为T样条曲面，进一步进行Volume Skinning
*/
void Test::test_generate_surfaces()
{
	string modelname = "isis";
	string filename = "../out/volume/" + modelname + "_bspline.txt";
	BsplineVolume volume;
	volume.readVolume(filename);
	

	int sample_num = 5;
	vector<NURBSSurface> nurbs(sample_num + 1);
	MatrixXd controlpoints; // 用于显示时设置相机位置和缩放大小
	// 采样生成B样条曲面
	for (int i = 0; i <= sample_num; i++) {
		double param = 1.0*i / sample_num;
		volume.get_isoparam_surface(nurbs[i], param, 'v');
		//nurbs[i].draw(Window::viewer, false, true);
		// 
		for (int j = 0; j < nurbs[i].controlPw.size(); j++) {
			int lastid = controlpoints.rows();
			controlpoints.conservativeResize(controlpoints.rows() + nurbs[i].controlPw[j].rows(), 3);
			for (int k = 0; k < nurbs[i].controlPw[j].rows(); k++) {
				controlpoints.row(lastid + k) = nurbs[i].controlPw[j].row(k);
			}
		}

	}
	
	vector<Mesh3d> tsplines(sample_num + 1);
	for (int i = 0; i <= sample_num; i++) {
		// venus_bspline.txt ---> 3e-3
		// tooth_bspline.txt ---> 3
		// isis_bspline.txt --> 5e-3
		TsplineSimplify(nurbs[i], tsplines[i], 20, 1e-2);
		cout << "number of nodes: " << tsplines[i].get_num() << endl;
		tsplines[i].setViewer(&Window::viewer);
		tsplines[i].draw(false, false, true, 0.01);
		Window::viewer.data_list[tsplines[i].id].set_colors(blue);
	}
	
	/*for (auto& surface : nurbs) {
		surface.draw(Window::viewer, false, true);
	}*/
	for (auto& data : Window::viewer.data_list) {
		data.set_face_based(true);
		data.show_lines = false;
		data.invert_normals = true;
	}
	cout << "controlpoints: " << controlpoints.rows() << endl;
	Window::viewer.core.align_camera_center(controlpoints);

	VolumeSkinning* method = new VolumePiaMethod(tsplines, 12, 1e-5);
	//VolumeSkinning* method = new VolumeSkinning(tsplines);
	method->setViewer(&Window::viewer);
	method->calculate();
	method->volume.saveVolume("../out/volume/" + modelname + "_skinning");
	method->volume.saveAsHex("../out/volume/" + modelname + "_skinning", 0.02);
	VolumeRender render(&method->volume, false, false, true, 0.01);
	begin = clock();
	render.launch();
	end = clock();
	cout << "time for drawing tspline volume: " << (end - begin) / CLOCKS_PER_SEC << "s" << endl;
	/*Window w;
	w.launch();*/
	
}

void Test::test_nurbs()
{
	// bezier curve
	/*NURBSCurve nurbs;
	nurbs.loadNURBS("../out/nurbs/bezier.cpt");
	nurbs.draw(Window::viewer);*/

	// bezier surface
	/*NURBSSurface nurbs;
	nurbs.loadNURBS("../out/nurbs/beziersurface.cpt");
	nurbs.draw(Window::viewer);*/

	// uniform bspline cruve
	NURBSCurve nurbs;
	nurbs.loadNURBS("../out/nurbs/bsplinecurve.cpt");
	nurbs.draw(Window::viewer);

	// uniform bspline surface 
	/*NURBSSurface nurbs;
	nurbs.loadNURBS("../out/nurbs/bsplinesurface.cpt");
	nurbs.draw(Window::viewer);*/

	// nurbs surface
	/*NURBSSurface nurbs;
	nurbs.loadNURBS("../out/nurbs/torus.cptw");
	nurbs.draw(Window::viewer);*/

	Window w;
	w.launch();

	
}

void Test::test_TsplineVolume() {
	TsplineVolume* volume = new TsplineVolume();
	volume->readVolume("../out/volume/isis_skinning.vol");
	
	VolumeRender render(volume, false, false, true, 0.01);
	begin = clock();
	render.launch();
	end = clock();
	cout << "time for drawing tspline volume: " << (end - begin) / CLOCKS_PER_SEC << "s" << endl;
}
void Test::test_BsplineVolume()
{
	BsplineVolume volume;
	begin = clock();
	volume.readVolume("../out/volume/tooth_bspline.txt");
	//volume.readVolume("../out/volume/balljoint_bspline.txt");
	//volume.readVolume("../out/volume/isis_bspline.txt");
	//volume.readVolume("../out/volume/moai_bspline.txt");
	//volume.readVolume("../out/volume/tooth_bspline.txt");
	VolumeRender render(&volume, false, false, true, 0.01);
	/*volume.saveAsHex("../out/volume/tooth_bspline", 0.1);
	volume.saveVolume("../out/volume/tooth_bspline");*/
	
	
	render.launch();
	end = clock();
	cout << "time passed: " << (end - begin) / CLOCKS_PER_SEC << "s" << endl;
}
void Test::test_Mesh() {

	//Mesh3d* mesh = new Mesh3d();
	//mesh->loadMesh("../out/tspline/simpleMesh2.cfg");
	//Mesh3d* meshcopy = new Mesh3d(*mesh); // deep copy
	//meshcopy->saveMesh("../out/tspline/simpleMesh2_copy");
	//delete mesh;
	//mesh = NULL;
	//MeshRender render(meshcopy);
	//render.launch();
	Mesh3d mesh;
	mesh.loadMesh("../out/tspline/venus_skinning.cfg");
	MeshRender render(&mesh, false, true, true);
	begin = clock();
	render.launch();
	end = clock();
	cout << "time for drawing tspline: " << (end - begin) / CLOCKS_PER_SEC << "s" << endl;

}

void Test::test_VolumeSkinning()
{
	vector<Mesh3d> surfaces(3);
	surfaces[0].loadMesh("../out/tspline/simpleMesh1.cfg");
	surfaces[1].loadMesh("../out/tspline/simpleMesh2.cfg");
	surfaces[2].loadMesh("../out/tspline/simpleMesh3.cfg");
	surfaces[0].setViewer(&Window::viewer);
	surfaces[1].setViewer(&Window::viewer);
	surfaces[2].setViewer(&Window::viewer);

	surfaces[0].draw(false, false, true);
	surfaces[1].draw(false, false, true);
	surfaces[2].draw(false, false, true);
	Window::viewer.data_list[1].set_colors(blue);
	Window::viewer.data_list[2].set_colors(blue);
	Window::viewer.data_list[3].set_colors(blue);
	/*Window w;
	w.launch();*/
	VolumeSkinning* method = new VolumePiaMethod(surfaces,15);
	method->setViewer(&Window::viewer);
	cout << "1" << endl;
	method->calculate();
	cout << "2" << endl;
	VolumeRender w(&method->volume, false, false, true,0.01);
	w.launch();
}

void Test::test_Skinning()
{

	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../out/nurbs/circle.cptw");
	nurbs[1].loadNURBS("../out/nurbs/circle1_1.cptw");
	nurbs[2].loadNURBS("../out/nurbs/circle2_1.cptw");
	nurbs[3].loadNURBS("../out/nurbs/circle3.cptw");

	nurbs[0].draw(Window::viewer, false);
	nurbs[1].draw(Window::viewer, false);
	nurbs[2].draw(Window::viewer, false);
	nurbs[3].draw(Window::viewer, false);
	
	//Skinning* method = new MinJaeMethod(nurbs, 20, 50);
	Skinning* method = new PiaMethod(nurbs, 100);
	//Skinning* method = new NasriMethod(nurbs);
	//Skinning* method = new OptMethod(nurbs);
	//Skinning* method = new PiaMinJaeMethod(nurbs, 100);
	//Skinning* method = new PiaNasriMethod(nurbs, 100);

	method->setViewer(&Window::viewer);
	method->calculate();
	Mesh3d* mesh = &(method->tspline);
	cout << "num of nodes: " << mesh->get_num() << endl;
	
	//mesh->saveMesh("../out/tspline/simpleMesh4");
	MeshRender render(mesh);
	render.launch();


}
void Test::test_DerOfNurbs() {
	NURBSCurve nurbs;
	nurbs.loadNURBS("../out/nurbs/circle.cptw");
	cout << "controlpw: \n" << nurbs.controlPw << endl;
	for (int i = 0; i <= 10; i++) {
		double u = 1.0*i / 10;
		MatrixXd point = MatrixXd::Zero(1, 3);
		point.row(0) = nurbs.eval(u);
		RowVector3d du = RowVector3d::Zero(); // 默认不是0

		for (int j = 0; j <= nurbs.n; j++) {
			double a = t_mesh::DersBasis(nurbs.knots, u, j, 3)(1);
			if (i == 0) {
				cout << "a: " << a << endl;
			}
			du += nurbs.controlPw.row(j) * a;
			if (i == 0) {
				cout << "du: \n" << du << endl;
			}
		}
		if (i == 0) {
			cout << "du: \n" << du << endl;
		}
		du.normalize();
		if (i == 0) {
			cout << "du: \n" << du << endl;
		}
		MatrixXd endpoint(1, 3);
		endpoint.row(0) = point.row(0) + du;
		Window::viewer.data().add_points(point, blue);
		Window::viewer.data().add_edges(point, endpoint, green);
	}
	nurbs.draw(Window::viewer);
	Window w;
	w.launch();
}
void Test::test_Lspia() {
	MatrixXd points;
	t_mesh::loadpoints("../out/points/helix.dat", points);
	
	NURBSCurve fit;
	fit.lspiafit(points, 10);

	fit.draw(Window::viewer, true,true,0.001);

	Window::viewer.data().add_points(points, red);

	Window w;
	w.launch();
}
void Test::test_Array() {
	t_mesh::Array<double, 5> A;
	A.input(cin);
	A = 2.0*A;
	A.output(cout);
	cout << endl;
	A = A * 2.0;
	A.output(cout);
	cout << endl;
}
void Test::test_Integral() {
	double a = 2.0;
	auto lambda = [a](double u, double v)-> double {
		return a*sin(u + v);
	};
	double res = OptMethod::integral(lambda);
	cout << "ingegral: " << res << endl;
	cout << "real: " << 2 * (2 * sin(1) - sin(2)) << endl;
}

void Test::test_Basis() {

	t_mesh::Array<double, 5> A;
	A.input(cin);
	A.output(cout);
	double t = 0.0;
	cin >> t;

	cout << "basis: " << t_mesh::Basis(A.toVectorXd(), t) << endl;
	cout << "basis1: " << Basis1(A.toVectorXd(), t) << endl;
}
void Test::test_Derivative() {
	Eigen::VectorXd knots(11);
	knots << 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5;
	double t = 2.5;
	cout << "derivative: \n" << DersBasis(knots, t, 4, 2) << endl;

	t_mesh::Array<double, 5> A;
	A.input(cin);
	A.output(cout);
	t = 0.0;
	cin >> t;

	cout << "derivative basis: \n" << t_mesh::DersBasis(A.toVectorXd(), t) << endl;

}
double Test::myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
	if (!grad.empty()) {
		grad[0] = 0.0;
		grad[1] = 0.5 / sqrt(x[1]);
	}
	return sqrt(x[1]);
}

double Test::myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
	my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
	double a = d->a, b = d->b;
	if (!grad.empty()) {
		grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
		grad[1] = -1.0;
	}
	return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}

void Test::test_Nlopt() {
	nlopt::opt opt(nlopt::LD_MMA, 2);
	std::vector<double> lb(2);
	lb[0] = -HUGE_VAL; lb[1] = 0;
	opt.set_lower_bounds(lb);
	opt.set_min_objective(myfunc, NULL);
	my_constraint_data data[2] = { { 2,0 },{ -1,1 } };
	opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
	opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
	opt.set_xtol_rel(1e-4);
	std::vector<double> x(2);
	x[0] = 1.234; x[1] = 5.678;
	double minf;

	try {
		nlopt::result result = opt.optimize(x, minf);
		std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
			<< std::setprecision(10) << minf << std::endl;
	}
	catch (std::exception &e) {
		std::cout << "nlopt failed: " << e.what() << std::endl;
	}
}