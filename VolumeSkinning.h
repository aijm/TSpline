#ifndef VOLUMESKINNING_H
#define VOLUMESKINNING_H

#include "TsplineVolume.h"
#include "NURBSCurve.h"
class VolumeSkinning {
public:
	VolumeSkinning(const vector<Mesh3d>& _surfaces) :surfaces(_surfaces) {
		surfaces_num = surfaces.size();
	}
	
	virtual void parameterize();    // 曲线参数化得到w方向参数 w_params
	virtual void init();        // 根据T样条曲面初始化体网格
	virtual void insert();      // 在体网格中插入中间面
	virtual void calculate();       // 计算流程

	// update coordinates of control points by the formula from (nasri 2012)
	// aX + bW + cY = V
	void update();

	// 用于skinning过程中显示中间结果
	void setViewer(Viewer* _viewer) {
		viewer = _viewer;
	}
private:
	Point3d centerOfmesh(const Mesh3d& mesh);
public:
	TsplineVolume volume;
protected:
	vector<Mesh3d> surfaces;
	int surfaces_num;
	Eigen::VectorXd w_params;
	Viewer* viewer;
};


#endif // !VOLUMESKINNING_H

