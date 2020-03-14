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
#include "PiaMinJaeMethod.h"
#include "PiaNasriMethod.h"
#include "TsplineVolume.h"
#include "BsplineVolume.h"
#include "NURBSSurface.h"
#include "VolumeSkinning.h"
#include "VolumePiaMethod.h"
#include "NurbsPia.h"

#include "MeshRender.h"
#include "VolumeRender.h"

class Test {
public:
	// surface skinning
	static void test_circle_skinning();
	static void test_venus_skinning();
	static void test_venus_skinning_helper_points();
	static void test_Bsurface_skinning();
	static void test_chess_skinning();
	static void test_ring_skinning();
	static void test_helicoidal_skinning();
	static void test_bonnet_skinning();
	static void test_door_skinning();
	static void test_face_skinning();

	// volume skinning
	static void test_sample_VolumeSkinning(string modelname, double simplifyEps, int sample_num = 5, char dir = 'v');
	static void test_VolumeSkinning(string modelname, double simplifyEps);


	// draw curve, surface, volume
	static void test_nurbs();
	static void test_TsplineVolume(const string& modelname, bool reverse = false);
	static void test_BsplineVolume(string modelname, double ratio = 0.01, bool reverse = false);
	static void test_Mesh();
	static void test_DrawMultiVolume();

	// helper function
	static void test_nurbscurve_interpolate_optimize();
	static void test_tspline_normal();
	static void test_getsurface_fromvolume();

	static void test_sample_fitbsplinesolid(string modelname, double simplifyEps);
	// 基于优化生成一个拟合给定数据点的有效B样条体
	static void test_fitbsplinesolid(string modelname, double simplifyEps);
	static void test_Nurbs_curvature();
	// B样条曲面简化为T样条曲面
	static void test_TsplineSimplify();

	// 读取另一种格式表示的nurbs曲面
	static void load_nurbs_surface(NURBSSurface& surface, string filename);
	// 保存另一种格式表示的nurbs曲面
	static void save_nurbs_surface(const NURBSSurface& surface, string filename);
	static void test_load_nurbs_surface();
	static void test_save_nurbs_surface();
	static void test_save_quadObj();

	static void test_DerOfNurbs();
	static void test_Lspia();
	static void test_Array();
	static void test_Basis();
	static void test_Derivative();

private:
	static clock_t begin;
	static clock_t end;

};

#endif // !TEST_H

