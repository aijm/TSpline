#ifndef UTILITY_H
#define UTILITY_H

#include<iostream>
#include<list>
#include<vector>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/gaussian_curvature.h>
#include <igl/jet.h>
#include <igl/parula.h>
#include <igl/principal_curvature.h>
#include "NURBSSurface.h"

namespace t_mesh{
    template<class T,int num>
    struct Array;

    template<class T>
    class Node;
    
    template<class T>
    class Mesh;
	
    typedef Array<float,3>  Point3f;
    typedef Array<double,3> Point3d;
    typedef Array<int,2>   Parameter;

	typedef Mesh<Point3d> Mesh3d;
	typedef Mesh<Point3d> Mesh3f;

	const Eigen::RowVector3d   red(1, 0, 0);
	const Eigen::RowVector3d green(0, 1, 0);
	const Eigen::RowVector3d  blue(0, 0, 1);
	const Eigen::RowVector3d  yellow(1, 1, 0);
	const Eigen::RowVector3d white(1, 1, 1);
	const Eigen::RowVector3d black(0, 0, 0);
	const Eigen::RowVector3d deeppink(255, 20, 147);

    // Blending function N[s0,s1,s2,s3,s4](p)
    template<class T>
    double B(const T& s,double p); 

	template<class T, int num>
	void array2matrixd(const Array<T, num> &a, Eigen::MatrixXd &m);

	int FindSpan(const Eigen::MatrixXd &knots, double t, int p = 3);

	// Blending function N[s0,s1,s2,s3,s4](p)
	double Basis1(const Eigen::MatrixXd &knotvector, double t, int i = 0, int p = 4);

	double Basis(const Eigen::MatrixXd &knots, double t, int i = 0, int p = 3);

	// Berivative of Blending function N[s0,s1,s2,s3,s4](t)
	Eigen::RowVectorXd DersBasis(const Eigen::MatrixXd &knots, double t, int i = 0, int p = 3);

	bool loadpoints(std::string name, Eigen::MatrixXd &mat);

	bool savepoints(std::string name, const Eigen::MatrixXd &mat);

	void vec_insert(Eigen::VectorXd &vec, double t);

	void TsplineSimplify(const NURBSSurface& surface, Mesh3d& tspline, int maxIterNum = 20, double eps = 1e-5);
	
};

namespace t_mesh{

    template<class T>
    double B(const T& s,double p){
        if(p<=s[0])
            return 0.0;
        if(p<=s[1])
            return 1.0*(p-s[0])*(p-s[0])*(p-s[0])/(s[1]-s[0])/(s[2]-s[0])/(s[3]-s[0]);
        else if(p<=s[2])
            return 1.0*(p-s[0])*(p-s[0])*(s[2]-p)/(s[2]-s[1])/(s[2]-s[0])/(s[3]-s[0])+
                1.0*(p-s[0])*(p-s[1])*(s[3]-p)/(s[2]-s[1])/(s[3]-s[0])/(s[3]-s[1])+
                1.0*(p-s[1])*(p-s[1])*(s[4]-p)/(s[2]-s[1])/(s[4]-s[1])/(s[3]-s[1])
                ;
        else if(p<=s[3])
            return 1.0*(p-s[0])*(s[3]-p)*(s[3]-p)/(s[3]-s[2])/(s[3]-s[1])/(s[3]-s[0])+
                1.0*(p-s[1])*(s[4]-p)*(s[3]-p)/(s[3]-s[2])/(s[4]-s[1])/(s[3]-s[1])+
                1.0*(p-s[2])*(s[4]-p)*(s[4]-p)/(s[3]-s[2])/(s[4]-s[2])/(s[4]-s[1])
                ;
        else if(p<s[4])
            return 1.0*(s[4]-p)*(s[4]-p)*(s[4]-p)/(s[4]-s[3])/(s[4]-s[2])/(s[4]-s[1]);
        return 0.0;
    }

	template<class T, int num>
	void array2matrixd(const Array<T, num> &a, Eigen::MatrixXd &m) {
		m.resize(1, num);
		for (int i = 0; i < num; i++) {
			m(0, i) = 1.0*a[i];
		}
	}

};



#include"array.hpp"
#include"node.hpp"
#include"mesh.hpp"

#endif
