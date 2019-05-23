//#include"utility.h"
#include<fstream>
#include<sstream>
#include<cstdlib>
#include "NasriMethod.h"
#include "MinJaeMethod.h"
#include "PiaMethod.h"
#include "window.h"

using namespace std;
using namespace Eigen;
using namespace window;

void testArray() {

	t_mesh::Array<double, 5> A;
	A.input(cin);
	A.output(cout);
	double t = 0.0;
	cin >> t;
	
	cout << "basis: " << t_mesh::Basis(A.toVectorXd(), t) << endl;

}
void testMesh(igl::opengl::glfw::Viewer &viewer) {
	mesh.loadMesh("../simpleMesh.cfg");
	mesh.draw(showmesh, showpolygon, showsurface);
}

void testSkinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../circle.cptw");
	nurbs[1].loadNURBS("../circle1_1.cptw");
	nurbs[2].loadNURBS("../circle2_1.cptw");
	nurbs[3].loadNURBS("../circle3.cptw");


	nurbs[0].draw(viewer,false);
	nurbs[1].draw(viewer,false);
	nurbs[2].draw(viewer,false);
	nurbs[3].draw(viewer, false);

	PiaMethod method(nurbs,1000);
	method.setViewer(&viewer);
	method.calculate();
	mesh = method.tspline;
	mesh.setViewer(&viewer);
	mesh.draw(false, true, true);

	cout << "num of nodes: " << mesh.get_num() << endl;
	//mesh.saveMesh("../simpleMesh1");
}


void testpia_skinning(igl::opengl::glfw::Viewer &viewer)
{
	vector<NURBSCurve> nurbs(4);
	nurbs[0].loadNURBS("../circle.cptw");
	nurbs[1].loadNURBS("../circle1.cptw");
	nurbs[2].loadNURBS("../circle2.cptw");
	nurbs[3].loadNURBS("../circle3.cptw");


	nurbs[0].draw(viewer, false);
	nurbs[1].draw(viewer, false);
	nurbs[2].draw(viewer, false);
	nurbs[3].draw(viewer, false);
	//viewer.core.align_camera_center(nurbs[1].controlPw);



	//cout << "knots: " << nurbs[0].knots.transpose() << endl;
	//mesh.pia_skinning(nurbs, viewer, 1000);

	mesh.draw(false, true, true);
	//mesh.saveMesh("simpleMesh");
}

int main(int argc,char** argv){
	//testArray();
	window::init();
	testSkinning(viewer);
	window::launch();
    return 0;
}




