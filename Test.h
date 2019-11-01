#ifndef TEST_H
#define TEST_H

#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
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
#include "BsplineVolume.h"
#include "NURBSSurface.h"
#include "VolumeSkinning.h"
#include "VolumePiaMethod.h"
#include "NurbsPia.h"
#include <nlopt.hpp>

#include "MeshRender.h"
#include "VolumeRender.h"

typedef struct {
	double a, b;
} my_constraint_data;


class Test {
public:
	// surface skinning
	static void test_circle_skinning();
	static void test_venus_skinning();
	static void test_venus_skinning_helper_points();
	static void test_Bsurface_skinning();
	static void test_chess_skinning();
	static void test_ring_skinning();

	// volume skinning
	static void test_sample_VolumeSkinning(string modelname, double simplifyEps, int sample_num = 5, char dir = 'v');
	static void test_VolumeSkinning(string modelname, double simplifyEps);


	// draw curve, surface, volume
	static void test_nurbs();
	static void test_TsplineVolume();
	static void test_BsplineVolume(string modelname, double ratio = 0.01, bool reverse = false);
	static void test_Mesh();

	// helper function
	static void test_nurbscurve_interpolate_optimize();
	static void test_tspline_normal();
	static void test_getsurface_fromvolume();

	static void test_fitbsplinesolid(string modelname, double simplifyEps);
	static void test_Nurbs_curvature();
	static void test_TsplineSimplify();

	static void load_nurbs_surface(NURBSSurface& surface, string filename);
	static void save_nurbs_surface(const NURBSSurface& surface, string filename);
	static void test_load_nurbs_surface();
	static void test_save_nurbs_surface();

	static void test_DerOfNurbs();
	static void test_Lspia();
	static void test_Array();
	static void test_Integral();
	static void test_Basis();
	static void test_Derivative();
	static void test_Nlopt();

private:
	static double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
	static double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data);
	
private:
	static clock_t begin;
	static clock_t end;

};

#endif // !TEST_H

