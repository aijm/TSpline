#include "window.h"
namespace window {
	using namespace std;

	t_mesh::Mesh3d mesh;
	double resolution = 0.01;
	bool showmesh = false;
	bool showpolygon = true;
	bool showsurface = true;
	// Customize the menu
	double doubleVariable = 0.1f; // Shared between two menus

	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	
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

				mesh.draw(showmesh, showpolygon, showsurface);
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
		else if (key == 'M') {
			viewer.data().clear();
			showmesh = !showmesh;
			mesh.draw(showmesh, showpolygon, showsurface);
		}
		else if (key == 'P') {
			viewer.data().clear();
			showpolygon = !showpolygon;
			mesh.draw(showmesh, showpolygon, showsurface);
		}
		else if (key == 'S') {
			viewer.data().clear();
			showsurface = !showsurface;
			mesh.draw(showmesh, showpolygon, showsurface);
		}
		return false;
	}

	void draw_viewer_menu()
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
	}

	void draw_custom_window()
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
	}

	void init()
	{
		viewer.plugins.push_back(&menu);
		// Add content to the default menu window
		menu.callback_draw_viewer_menu = &draw_viewer_menu;

		// Draw additional windows
		menu.callback_draw_custom_window = &draw_custom_window;

		// key binding
		viewer.callback_key_down = &key_down;

		mesh.setViewer(&viewer);
	}

	void launch()
	{
		viewer.launch();
	}

};