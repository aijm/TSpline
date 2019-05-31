#ifndef VOLUMERENDER_H
#define VOLUMERENDER_H

#include "Window.h"
#include "Volume.h"
class VolumeRender : public Window {

public:
	VolumeRender(Volume* _volume,
		bool _showmesh = false, bool _showpolygon = true,
		bool _showvolume = true, double _resolution = 0.1)
		:volume(_volume), showmesh(_showmesh), showpolygon(_showpolygon),
		showvolume(_showvolume), resolution(_resolution)
	{
		volume->setViewer(&viewer);
	}

protected:
	void draw() override;

private:
	Volume* volume;
	bool showmesh;
	bool showpolygon;
	bool showvolume;
	double resolution;
};

#endif // !VOLUMERENDER_H

