#include"utility.h"
#include"draw.h"
#include<fstream>
#include<sstream>
//#include<boost/timer.hpp>
#include<cstdlib>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

using namespace std;
t_mesh::Mesh3d mesh;
double resolution = 0.01;
bool showmesh = false;
bool showpolygon = true;
bool showsurface = true;

void insert_loop(igl::opengl::glfw::Viewer &viewer) {
	double s = 0.0;
	double t = 0.0;
	while (true) {
		cout << "insert kont, format: s t" << endl;
		
		if (!(cin >> s >> t)) {
			cin.clear(); //clear the buffer
			cin.get();
			cout << "error! please use right format!" << endl;
			continue;
		}
		else {
			mesh.insert(s, t);
			//mesh.drawTmesh(viewer);
			//mesh.drawControlpolygon(viewer);
			//mesh.drawSurface(viewer);
			mesh.draw(viewer, showmesh, showpolygon, showsurface);
			return;
		}
	}
}
// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	if (key == 32) {
		viewer.data().clear();
		insert_loop(viewer);
	}
	else if(key == 'M'){
		viewer.data().clear();
		showmesh = !showmesh;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	else if (key == 'P') {
		viewer.data().clear();
		showpolygon = !showpolygon;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	else if (key == 'S') {
		viewer.data().clear();
		showsurface = !showsurface;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	return false;
}

void testArray() {
	/*t_mesh::Point3d A;
	A.input(cin);
	A.output(cout);*/
	/*Eigen::MatrixXd m;
	t_mesh::array2matrixd(A, m);
	cout << "tomatrix: \n" << m << endl;*/
	t_mesh::Array<double, 5> A;
	A.input(cin);
	A.output(cout);
	double t = 0.0;
	cin >> t;
	
	cout << "basis: " << t_mesh::Basis(A.toVectorXd(), t) << endl;

}
void testMesh(igl::opengl::glfw::Viewer &viewer) {
	
	
	//mesh.insert(1.5, 1.5);
	//mesh.insert(2.0, 1.5);
	////mesh.drawControlpolygon(viewer);
	////mesh.drawSurface(viewer, 0.01);
	//mesh.drawTmesh(viewer);
}

void testSkinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../circle.cptw");
	nurbs[1].loadNURBS("../circle1_1.cptw");
	nurbs[2].loadNURBS("../circle2_1.cptw");
	nurbs[3].loadNURBS("../circle3.cptw");


	nurbs[0].draw(viewer,false);
	nurbs[1].draw(viewer,false);
	nurbs[2].draw(viewer,false);
	nurbs[3].draw(viewer, false);
	//viewer.core.align_camera_center(nurbs[1].controlPw);



	//cout << "knots: " << nurbs[0].knots.transpose() << endl;
	mesh.skinning(nurbs, viewer);

	mesh.draw(viewer, false,true, true);
	//mesh.saveMesh("simpleMesh");
}


void testpia_skinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../circle.cptw");
	nurbs[1].loadNURBS("../circle1_1.cptw");
	nurbs[2].loadNURBS("../circle2_1.cptw");
	nurbs[3].loadNURBS("../circle3.cptw");


	nurbs[0].draw(viewer, false);
	nurbs[1].draw(viewer, false);
	nurbs[2].draw(viewer, false);
	nurbs[3].draw(viewer, false);
	//viewer.core.align_camera_center(nurbs[1].controlPw);



	//cout << "knots: " << nurbs[0].knots.transpose() << endl;
	mesh.pia_skinning(nurbs, viewer, 200);

	mesh.draw(viewer, false, true, true);
	//mesh.saveMesh("simpleMesh");
}

int main(int argc,char** argv){
	igl::opengl::glfw::Viewer viewer;
	//testArray();
	testpia_skinning(viewer);
	//testSkinning(viewer);
	/*if (mesh.loadMesh("../surface1.cfg") == 0) {
		cout << "num of nodes: " << mesh.get_num() << endl;
	}
	else {
		cout << "load failing!" << endl;
	}*/
	
	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Customize the menu
	double doubleVariable = 0.1f; // Shared between two menus

								  // Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::InputDouble("double", &doubleVariable, 0, 0, "%.4f");

			// ... or using a custom callback
			static bool boolVariable = true;
			if (ImGui::Checkbox("bool", &boolVariable))
			{
				// do something
				std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
			}

			// Expose an enumeration type
			enum Orientation { Up = 0, Down, Left, Right };
			static Orientation dir = Up;
			ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

			// We can also use a std::vector<std::string> defined dynamically
			static int num_choices = 3;
			static std::vector<std::string> choices;
			static int idx_choice = 0;
			if (ImGui::InputInt("Num letters", &num_choices))
			{
				num_choices = std::max(1, std::min(26, num_choices));
			}
			if (num_choices != (int)choices.size())
			{
				choices.resize(num_choices);
				for (int i = 0; i < num_choices; ++i)
					choices[i] = std::string(1, 'A' + i);
				if (idx_choice >= num_choices)
					idx_choice = num_choices - 1;
			}
			ImGui::Combo("Letter", &idx_choice, choices);

			// Add a button
			if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
			{
				std::cout << "Hello\n";
			}
		}
	};

	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"New Window", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);


		// Expose the same variable directly ...
		ImGui::PushItemWidth(-80);
		ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
		ImGui::PopItemWidth();

		static std::string str = "bunny";
		ImGui::InputText("Name", str);

		ImGui::End();
	};

	viewer.callback_key_down = &key_down;
	//mesh.draw(viewer, showmesh, showpolygon, showsurface);
	viewer.launch();
	
    return 0;
}

