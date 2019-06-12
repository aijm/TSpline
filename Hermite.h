#ifndef HERMITE_H
#define HERMITE_H
#include<iostream>
#include<random>
#include<vector>
#include <igl/opengl/glfw/Viewer.h>
using namespace std;
using namespace igl::opengl::glfw;
using namespace Eigen;
struct Hermite
{
	//RowVectorXd eval(double t);
	VectorXd x;
	MatrixXd y;
	MatrixXd dy;
	double d2_y0 = 0;
	double d2_yn = 0;
	
};
#endif // !HERMITE_H

