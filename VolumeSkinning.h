#pragma once
#include "TsplineVolume.h"
class VolumeSkinning {
public:
	VolumeSkinning(const vector<Mesh3d>& _surfaces):surfaces(_surfaces){}


	TsplineVolume volume;
protected:

	vector<Mesh3d> surfaces;
};
