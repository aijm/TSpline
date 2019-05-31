#include "VolumeRender.h"

void VolumeRender::draw()
{
	volume->draw(showmesh, showpolygon, showvolume, resolution);
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
			volume->draw(showmesh, showpolygon, showvolume, resolution);
		}
		else if (key == 'P') {
			viewer.data().clear();
			showpolygon = !showpolygon;
			volume->draw(showmesh, showpolygon, showvolume, resolution);
		}
		else if (key == 'S') {
			viewer.data().clear();
			showvolume = !showvolume;
			volume->draw(showmesh, showpolygon, showvolume, resolution);
		}
		return false;
	};
}
