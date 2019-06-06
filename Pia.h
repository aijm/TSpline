#include "FitPoint.hpp"
#include "utility.h"
using namespace std;
using namespace t_mesh;
class Pia{
public:
	Pia(std::vector<FitPoint2D> _fitPoints, int _maxIterNum = 100, double _eps = 1e-5)
		:fitPoints(_fitPoints),maxIterNum(_maxIterNum), eps(_eps) {

	}
	void init();		// 根据NUUBSCurve初始化T-preimage
	void adapative_insert();		// 按一定规则在误差大的地方插入节点，局部加细
	void calculate();           // 计算流程

	void fit();
	void pia();
public:
	Mesh3d tspline;
private:
	const int maxIterNum;
	const double eps;
	double error;
	vector<FitPoint2D> fitPoints;
};