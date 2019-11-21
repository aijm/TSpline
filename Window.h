#ifndef WINDOW_H
#define WINDOW_H

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

class Window {
public:
	void launch();
	~Window();

protected:
	virtual void init();
	virtual void draw();
public:
	static igl::opengl::glfw::Viewer viewer;

private:
	static double doubleVariable; // Shared between two menus
	static igl::opengl::glfw::imgui::ImGuiMenu menu;
};

#endif // !WINDOW_H

