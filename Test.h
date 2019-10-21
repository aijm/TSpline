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
	static void test_nurbs_pia();
	static void load_nurbs_surface(NURBSSurface& surface, string filename);
	static void save_nurbs_surface(const NURBSSurface& surface, string filename);
	static void test_load_nurbs_surface();
	static void test_save_nurbs_surface();
	static void test_chess();
	static void test_generate_surfaces();
	static void test_Nurbs_curvature();
	static void test_TsplineSimplify();
	static void test_generate_curves();
	static void test_generate_curves1();
	static void test_generate_curves2();
	static void test_nurbs();
	static void test_TsplineVolume();
	static void test_BsplineVolume();
	static void test_Mesh();
	static void test_VolumeSkinning();
	static void test_Skinning();
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

