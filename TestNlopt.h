#pragma once
#include <iomanip>
#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include "OptMethod.h"

typedef struct {
	double a, b;
} my_constraint_data;

class TestNlopt {
public:
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
		cout << "real: " << 2* (2 * sin(1) - sin(2)) << endl;
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