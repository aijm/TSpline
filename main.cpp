#include<fstream>
#include<sstream>
#include<cstdlib>
#include "NasriMethod.h"
#include "MinJaeMethod.h"
#include "PiaMethod.h"
#include "OptMethod.h"
#include "PiaMinJaeMethod.h"
#include "PiaNasriMethod.h"
#include "window.h"

#include "TestNlopt.h"

using namespace std;
using namespace Eigen;
using namespace window;



void testMesh(igl::opengl::glfw::Viewer &viewer) {
	mesh.loadMesh("../simpleMesh.cfg");
	mesh.draw(showmesh, showpolygon, showsurface);
}

void testSkinning()
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
	Skinning* method = new MinJaeMethod(nurbs, 100, 10);
	//Skinning* method = new PiaMethod(nurbs, 1000);
	//Skinning* method = new NasriMethod(nurbs);
	//Skinning* method = new OptMethod(nurbs);
	//Skinning* method = new PiaMinJaeMethod(nurbs, 1000);
	//Skinning* method = new PiaNasriMethod(nurbs, 1000);

	method->setViewer(&viewer);
	method->calculate();
	mesh = method->tspline;
	mesh.setViewer(&viewer);
	mesh.draw(false, true, true);

	cout << "num of nodes: " << mesh.get_num() << endl;
	//mesh.saveMesh("../simpleMesh1");
}
void test_derOfNurbs() {
	NURBSCurve nurbs;
	nurbs.loadNURBS("../circle.cptw");
	cout << "controlpw: \n" << nurbs.controlPw << endl;
	for (int i = 0; i <= 10; i++) {
		double u = 1.0*i / 10;
		MatrixXd point = MatrixXd::Zero(1, 3);
		point.row(0) = nurbs.eval(u);
		RowVector3d du = RowVector3d::Zero(); // Ä¬ÈÏ²»ÊÇ0

		for (int j = 0; j <= nurbs.n; j++) {
			double a = t_mesh::DersBasis(nurbs.knots, u, j, 3)(1);
			if (i == 0) {
				cout << "a: " << a << endl;
			}
			du += nurbs.controlPw.row(j) * a;
			if (i == 0) {
				cout << "du: \n" << du << endl;
			}
		}
		if (i == 0) {
			cout << "du: \n" << du << endl;
		}
		du.normalize();
		if (i == 0) {
			cout << "du: \n" << du << endl;
		}
		MatrixXd endpoint(1, 3);
		endpoint.row(0) = point.row(0) + du;
		viewer.data().add_points(point, blue);
		viewer.data().add_edges(point, endpoint, green);
	}
	nurbs.draw(viewer);
}
void test_lspia() {
	NURBSCurve nurbs;
	nurbs.loadNURBS("../circle.cptw");
	const int sampleNum = 100;
	MatrixXd points(sampleNum + 1, nurbs.controlPw.cols());
	VectorXd params(points.rows());
	for (int i = 0; i <= sampleNum; i++) {
		params(i) = 1.0*i / sampleNum;
		points.row(i) = nurbs.eval(params(i));
	}
	NURBSCurve fit;
	fit.lspiafit(points, params,nurbs.controlPw.rows(), nurbs.knots, 1000);

	nurbs.draw(viewer);
	fit.draw(viewer);
}




int main(int argc,char** argv){
	/*TestNlopt::test_Array();
	TestNlopt::test_nlopt();
	TestNlopt::test_Integral();*/

	//TestNlopt::testBasis();
	//TestNlopt::testDerivative();
	
	window::init();
	//test_derOfNurbs();
	//test_lspia();
	testSkinning();
	window::launch();
    return 0;
}




