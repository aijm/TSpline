#ifndef TEST_H
#define TEST_H

#include <iomanip>
#include <iostream>
#include <vector>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include "NasriMethod.h"
#include "MinJaeMethod.h"
#include "PiaMethod.h"
#include "OptMethod.h"
#include "PiaMinJaeMethod.h"
#include "PiaNasriMethod.h"
#include "TsplineVolume.h"

#include <nlopt.hpp>

#include "window.h"

using namespace window;
typedef struct {
	double a, b;
} my_constraint_data;

class Test {
public:
	static void testTsplineVolume() {
		TsplineVolume volume;
	}
	static void testMesh() {
		mesh.loadMesh("../simpleMesh2.cfg");
		mesh.draw(showmesh, showpolygon, showsurface);
	}

	static void testSkinning()
	{
		vector<NURBSCurve> nurbs(4);
		nurbs[0].loadNURBS("../circle.cptw");
		nurbs[1].loadNURBS("../circle1.cptw");
		nurbs[2].loadNURBS("../circle2.cptw");
		nurbs[3].loadNURBS("../circle3.cptw");


		nurbs[0].draw(viewer, false);
		nurbs[1].draw(viewer, false);
		nurbs[2].draw(viewer, false);
		nurbs[3].draw(viewer, false);
		//Skinning* method = new MinJaeMethod(nurbs, 30, 0);
		//Skinning* method = new PiaMethod(nurbs, 1000);
		Skinning* method = new NasriMethod(nurbs);
		//Skinning* method = new OptMethod(nurbs);
		//Skinning* method = new PiaMinJaeMethod(nurbs, 1000);
		//Skinning* method = new PiaNasriMethod(nurbs, 1000);

		method->setViewer(&viewer);
		method->calculate();
		mesh = method->tspline;
		mesh.setViewer(&viewer);
		mesh.draw(false, true, true);

		cout << "num of nodes: " << mesh.get_num() << endl;
		for (auto node : mesh.nodes) {
			Point3d temp = node->data;
			temp[0] += 1;
			node->data = temp;
		}
		mesh.saveMesh("../simpleMesh2");
	}
	static void test_derOfNurbs() {
		NURBSCurve nurbs;
		nurbs.loadNURBS("../circle.cptw");
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
			viewer.data().add_points(point, blue);
			viewer.data().add_edges(point, endpoint, green);
		}
		nurbs.draw(viewer);
	}
	static void test_lspia() {
		NURBSCurve nurbs;
		nurbs.loadNURBS("../circle.cptw");
		const int sampleNum = 100;
		MatrixXd points(sampleNum + 1, nurbs.controlPw.cols());
		VectorXd params(points.rows());
		for (int i = 0; i <= sampleNum; i++) {
			params(i) = 1.0*i / sampleNum;
			points.row(i) = nurbs.eval(params(i));
		}
		NURBSCurve fit;
		fit.lspiafit(points, params, nurbs.controlPw.rows(), nurbs.knots, 1000);

		nurbs.draw(viewer);
		fit.draw(viewer);
	}
	static void test_Array() {
		t_mesh::Array<double, 5> A;
		A.input(cin);
		A = 2.0*A;
		A.output(cout);
		cout << endl;
		A = A * 2.0;
		A.output(cout);
		cout << endl;
	}
	static void test_Integral() {
		double a = 2.0;
		auto lambda = [a](double u, double v)-> double {
			return a*sin(u + v);
		};
		double res = OptMethod::integral(lambda);
		cout << "ingegral: " << res << endl;
		cout << "real: " << 2 * (2 * sin(1) - sin(2)) << endl;
	}

	static void testBasis() {

		t_mesh::Array<double, 5> A;
		A.input(cin);
		A.output(cout);
		double t = 0.0;
		cin >> t;

		cout << "basis: " << t_mesh::Basis(A.toVectorXd(), t) << endl;
		cout << "basis1: " << Basis1(A.toVectorXd(), t) << endl;
	}
	static void testDerivative() {
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
	static double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
	{
		if (!grad.empty()) {
			grad[0] = 0.0;
			grad[1] = 0.5 / sqrt(x[1]);
		}
		return sqrt(x[1]);
	}

	static double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
	{
		my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
		double a = d->a, b = d->b;
		if (!grad.empty()) {
			grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
			grad[1] = -1.0;
		}
		return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
	}

	static void test_nlopt() {
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

};

#endif // !TEST_H

