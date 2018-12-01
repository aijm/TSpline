#include"utility.h"
#include"draw.h"
#include<fstream>
#include<sstream>
//#include<boost/timer.hpp>
#include<cstdlib>

using namespace std;

void testArray() {
	t_mesh::Point3d A;
	A.input(cin);
	A.output(cout);
	Eigen::MatrixXd m;
	t_mesh::array2matrixd(A, m);
	cout << "tomatrix: \n" << m << endl;
}
void testMesh(igl::opengl::glfw::Viewer &viewer) {
	t_mesh::Mesh3d mesh;
	if (mesh.loadMesh("../surface1.cfg") == 0) {
		cout << "num of nodes: " << mesh.get_num() << endl;
	}
	else {
		cout << "load failing!" << endl;
	}
	//mesh.drawTmesh(viewer);
	mesh.drawControlpolygon(viewer);
	mesh.drawSurface(viewer, 0.1);
}
int main(int argc,char** argv){
	igl::opengl::glfw::Viewer viewer;
	//testArray();
	testMesh(viewer);
	viewer.launch();
    return 0;
}

