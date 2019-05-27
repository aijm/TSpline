#include<fstream>
#include<sstream>
#include<cstdlib>
#include "NasriMethod.h"
#include "MinJaeMethod.h"
#include "PiaMethod.h"
#include "window.h"

#include "TestNlopt.h"

using namespace std;
using namespace Eigen;
using namespace window;

void testBasis() {

	t_mesh::Array<double, 5> A;
	A.input(cin);
	A.output(cout);
	double t = 0.0;
	cin >> t;
	
	cout << "basis: " << t_mesh::Basis(A.toVectorXd(), t) << endl;
	cout << "basis1: " << Basis1(A.toVectorXd(), t) << endl;


}
void testDerivative() {
	Eigen::VectorXd knots(11);
	knots << 0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5;
	double t = 2.5;
	cout << "derivative: \n" << DersBasis(knots, t, 4, 2) << endl;
}
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
	Skinning* method = new MinJaeMethod(nurbs, 100, 0);
	//Skinning* method = new PiaMethod(nurbs, 1000);
	method->setViewer(&viewer);
	method->calculate();
	mesh = method->tspline;
	mesh.setViewer(&viewer);
	mesh.draw(false, true, true);

	cout << "num of nodes: " << mesh.get_num() << endl;
	//mesh.saveMesh("../simpleMesh1");
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
	//TestNlopt::test_nlopt();
	//testBasis();
	testDerivative();
	window::init();
	//test_lspia();
	testSkinning();
	window::launch();
    return 0;
}




