#ifndef UTILITY_H
#define UTILITY_H

#include<iostream>
#include<list>
#include<vector>
//#include<cv.h>
//#include<highgui.h>
#include <igl/opengl/glfw/Viewer.h>

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

    // Blending function N[s0,s1,s2,s3,s4](p)
    template<class T>
    double B(const T& s,double p); 
};

namespace t_mesh{
    // Blending function N[s0,s1,s2,s3,s4](p)
	double Basis(const Eigen::MatrixXd &knotvector, double t,int i=0, int p=4)
	{
		//cout << "knotvector:\n" << knotvector.transpose() << endl;
		//int p = knotvector.size() - 1;
		assert(p >= 1);
		if (p == 1) {
			if (t >= knotvector(i) && t < knotvector(i + 1)) {
				return 1.0;
			}
			else {
				return 0.0;
			}
		}

		double a = knotvector(i + p - 1) - knotvector(i);
		double b = knotvector(i + p) - knotvector(i + 1);
		a = (a == 0.0) ? 0.0 : (t - knotvector(i)) / a;
		b = (b == 0.0) ? 0.0 : (knotvector(i + p) - t) / b;
		return a*Basis(knotvector, t,i, p - 1) + b*Basis(knotvector, t, i + 1, p - 1);
	}




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

};

#include"array.hpp"
#include"node.hpp"
#include"mesh.hpp"

#endif
