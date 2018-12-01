#ifndef DRAW_H
#define DRAW_H
#include <igl/opengl/glfw/Viewer.h>

namespace t_mesh {
	const Eigen::RowVector3d   red(1, 0, 0);
	const Eigen::RowVector3d green(0, 1, 0);
	const Eigen::RowVector3d  blue(0, 0, 1);
	const Eigen::RowVector3d white(1, 1, 1);
	const Eigen::RowVector3d black(0, 0, 0);
	

	template<class T,int num>
	void array2matrixd(const Array<T, num> &a, Eigen::MatrixXd &m) {
		m.resize(1, num);
		for (int i = 0; i < num; i++) {
			m(0, i) = 1.0*a[i];
		}
	}
	

};

#endif // !DRAW_H

