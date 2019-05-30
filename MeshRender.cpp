#include "MeshRender.h"

void MeshRender::draw()
{
	mesh->draw(showmesh, showpolygon, showsurface, resolution);
	viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) ->bool
	{
		//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
		if (key == 32) {
			/*viewer.data().clear();
			insert_loop(viewer);*/
		}
		else if (key == 'M') {
			viewer.data().clear();
			showmesh = !showmesh;
			mesh->draw(showmesh, showpolygon, showsurface);
		}
		else if (key == 'P') {
			viewer.data().clear();
			showpolygon = !showpolygon;
			mesh->draw(showmesh, showpolygon, showsurface);
		}
		else if (key == 'S') {
			viewer.data().clear();
			showsurface = !showsurface;
			mesh->draw(showmesh, showpolygon, showsurface);
		}
		return false;
	};
}
