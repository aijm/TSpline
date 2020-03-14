#include "Test.h"

int main(int argc,char** argv){
	// surface skinning
	//Test::test_circle_skinning();
	//Test::test_venus_skinning();
	//Test::test_venus_skinning_helper_points();
	//Test::test_Bsurface_skinning();
	//Test::test_chess_skinning();
	//Test::test_ring_skinning();
	//Test::test_helicoidal_skinning();
	//Test::test_bonnet_skinning();
	//Test::test_door_skinning();
	//Test::test_face_skinning();

	// volume skinning
	//Test::test_sample_VolumeSkinning("tooth", 3);
	//Test::test_sample_VolumeSkinning("venus", 3e-3);
	//Test::test_sample_VolumeSkinning("isis", 5e-3);
	//Test::test_sample_VolumeSkinning("isis", 1e-2);
	//Test::test_sample_VolumeSkinning("moai", 0.04, 5, 'u');
	//Test::test_sample_VolumeSkinning("head", 0.01);
	//Test::test_VolumeSkinning("tooth", 3);
	//Test::test_VolumeSkinning("venus", 3e-3);
	//Test::test_VolumeSkinning("isis", 5e-3);
	//Test::test_VolumeSkinning("isis", 1e-2);
	//Test::test_VolumeSkinning("moai", 0.04);
	//Test::test_VolumeSkinning("moai_new", 0.03);
	//Test::test_VolumeSkinning("moai_fitbspline", 0.04);
	//Test::test_VolumeSkinning("Ssolid", 0.05);
	//Test::test_VolumeSkinning("head", 0.01);
	//Test::test_VolumeSkinning("duck", 2);
	//Test::test_VolumeSkinning("duck_new", 2);


	// draw curve, surface, volume
	//Test::test_nurbs();
	//Test::test_TsplineVolume("duck_new", true);
	//Test::test_TsplineVolume("tooth", false);
	//Test::test_TsplineVolume("Ssolid", true);
	//Test::test_TsplineVolume("moai_new", true);
	Test::test_BsplineVolume("tooth", 0.02, true);
	//Test::test_BsplineVolume("venus", 0.01, true);
	//Test::test_BsplineVolume("moai", 0.02, true);
	//Test::test_BsplineVolume("balljoint", 0.01, true);
	//Test::test_BsplineVolume("isis", 0.01, true);
	//Test::test_BsplineVolume("head", 0.02, true);
	//Test::test_BsplineVolume("head", 0.02, true);
	//Test::test_BsplineVolume("multiVolume_4", 0.02, true);
	//Test::test_BsplineVolume("moai_new", 0.02, true);
	//Test::test_BsplineVolume("duck", 0.02, true);
	//Test::test_BsplineVolume("duck_new", 0.02, true);
	//Test::test_Mesh();
	//Test::test_DrawMultiVolume();

	// helper function
	//Test::test_nurbscurve_interpolate_optimize();
	//Test::test_tspline_normal();
	//Test::test_getsurface_fromvolume();

	//Test::test_fitbsplinesolid("Ssolid", 0.05);
	//Test::test_fitbsplinesolid("moai_fitbspline", 0.04);
	//Test::test_sample_fitbsplinesolid("moai", 0.04);
	//Test::test_Nurbs_curvature();
	//Test::test_TsplineSimplify();

	
	//Test::test_load_nurbs_surface();
	//Test::test_save_nurbs_surface();
	//Test::test_save_quadObj();
	//Test::test_DerOfNurbs();
	//Test::test_Lspia();
	//Test::test_Array();
	//Test::test_Basis();
	//Test::test_Derivative();
    return 0;
}




