#pragma once
#include "utility.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

namespace window {
	extern t_mesh::Mesh3d mesh;
	extern double resolution;
	extern bool showmesh;
	extern bool showpolygon;
	extern bool showsurface;
	extern double doubleVariable;

	extern igl::opengl::glfw::Viewer viewer;
	// Attach a menu plugin
	extern igl::opengl::glfw::imgui::ImGuiMenu menu;

	void insert_loop(igl::opengl::glfw::Viewer &viewer);

	// callback_key_down
	bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);

	// callback_draw_viewer_menu
	void draw_viewer_menu();

	// callback_draw_custom_window
	void draw_custom_window();

	void init();
	void launch();

};
