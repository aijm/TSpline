#include "Test.h"

clock_t Test::begin;
clock_t Test::end;

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
	volume->readVolume("../out/test.vol");

	TsplineVolume volume_copy(*volume); // deep copy
	delete volume;
	volume = NULL;
	VolumeRender render(&volume_copy, false, false, true);
	volume_copy.saveVolume("../out/test_copy");
	render.launch();
}
void Test::test_BsplineVolume()
{
	BsplineVolume volume;
	begin = clock();
	volume.readVolume("../out/venus_bspline.txt");
	//volume.readVolume("../out/balljoint_bspline.txt");
	//volume.readVolume("../out/isis_bspline.txt");
	//volume.readVolume("../out/moai_bspline.txt");
	//volume.readVolume("../out/tooth_bspline.txt");
	VolumeRender render(&volume, false, false, true, 0.01);
	/*volume.saveAsHex("../out/tooth_bspline", 0.1);
	volume.saveVolume("../out/tooth_bspline");*/
	
	
	render.launch();
	end = clock();
	cout << "time passed: " << (end - begin) / CLOCKS_PER_SEC << "s" << endl;
}
void Test::test_Mesh() {

	Mesh3d* mesh = new Mesh3d();
	mesh->loadMesh("../out/simpleMesh2.cfg");
	Mesh3d* meshcopy = new Mesh3d(*mesh); // deep copy
	meshcopy->saveMesh("../out/simpleMesh2_copy.cfg");
	delete mesh;
	mesh = NULL;
	MeshRender render(meshcopy);
	render.launch();

}

void Test::test_VolumeSkinning()
{
	vector<Mesh3d> surfaces;
	//vector<Mesh3d> a(surfaces);
	VolumeSkinning vs(surfaces);

}

void Test::test_Skinning()
{

	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../out/nurbs/circle.cptw");
	nurbs[1].loadNURBS("../out/nurbs/circle1.cptw");
	nurbs[2].loadNURBS("../out/nurbs/circle2.cptw");
	nurbs[3].loadNURBS("../out/nurbs/circle3.cptw");


	nurbs[0].draw(Window::viewer, false);
	nurbs[1].draw(Window::viewer, false);
	nurbs[2].draw(Window::viewer, false);
	nurbs[3].draw(Window::viewer, false);
	//Skinning* method = new MinJaeMethod(nurbs, 100, 10);
	//Skinning* method = new PiaMethod(nurbs, 1000);
	Skinning* method = new NasriMethod(nurbs);
	//Skinning* method = new OptMethod(nurbs);
	//Skinning* method = new PiaMinJaeMethod(nurbs, 1000);
	//Skinning* method = new PiaNasriMethod(nurbs, 1000);

	method->setViewer(&Window::viewer);
	method->calculate();
	Mesh3d* mesh = &(method->tspline);
	cout << "num of nodes: " << mesh->get_num() << endl;
	for (auto node : mesh->nodes) {
		Point3d temp = node->data;
		temp[0] += 6;
		node->data = temp;
	}
	mesh->saveMesh("../out/simpleMesh2");
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
		RowVector3d du = RowVector3d::Zero(); // Ä¬ÈÏ²»ÊÇ0

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
	NURBSCurve nurbs;
	nurbs.loadNURBS("../out/nurbs/circle.cptw");
	const int sampleNum = 100;
	MatrixXd points(sampleNum + 1, nurbs.controlPw.cols());
	VectorXd params(points.rows());
	for (int i = 0; i <= sampleNum; i++) {
		params(i) = 1.0*i / sampleNum;
		points.row(i) = nurbs.eval(params(i));
	}
	NURBSCurve fit;
	fit.lspiafit(points, params, nurbs.controlPw.rows(), nurbs.knots, 1000);

	nurbs.draw(Window::viewer);
	fit.draw(Window::viewer);
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