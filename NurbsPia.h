#ifndef NURBSPIA_H
#define NURBSPIA_H

#include "NURBSSurface.h"

#include "PiaMethod.h"
class NurbsPia {
public:
	NurbsPia::NurbsPia(const NURBSSurface & _surface, int _maxIterNum = 20, double _eps = 1e-5, int _num = 20)
		: surface(_surface), maxIterNum(_maxIterNum), eps(_eps), num(_num)
	{
		init();
	}
	void init();
	void fit();
	void pia();
	NURBSSurface calculate();
	
public:
	const int num;
	const int maxIterNum;
	const double eps;
	double error;
	vector<MatrixXd> eval;
	vector<MatrixXd> points;
	NURBSSurface surface;
	NURBSSurface pia_surface;
};

#endif // !NURBSPIA_H

