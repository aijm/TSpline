#include"utility.h"
#include"draw.h"
#include<fstream>
#include<sstream>
//#include<boost/timer.hpp>
#include<cstdlib>


using namespace std;
t_mesh::Mesh3d mesh;
double resolution = 0.01;
bool showmesh = true;
bool showpolygon = true;
bool showsurface = true;

void insert_loop(igl::opengl::glfw::Viewer &viewer) {
	double s = 0.0;
	double t = 0.0;
	while (true) {
		cout << "insert kont, format: s t" << endl;
		
		if (!(cin >> s >> t)) {
			cin.clear(); //锟斤拷栈锟斤拷锟斤拷锟�
			cin.get();
			cout << "error! please use right format!" << endl;
			continue;
		}
		else {
			mesh.insert(s, t);
			//mesh.drawTmesh(viewer);
			//mesh.drawControlpolygon(viewer);
			//mesh.drawSurface(viewer);
			mesh.draw(viewer, showmesh, showpolygon, showsurface);
			return;
		}
	}
}
// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
	//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	if (key == 32) {
		viewer.data().clear();
		insert_loop(viewer);
	}
	else if(key == 'M'){
		viewer.data().clear();
		showmesh = !showmesh;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	else if (key == 'P') {
		viewer.data().clear();
		showpolygon = !showpolygon;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	else if (key == 'S') {
		viewer.data().clear();
		showsurface = !showsurface;
		mesh.draw(viewer, showmesh, showpolygon, showsurface);
	}
	return false;
}

void testArray() {
	t_mesh::Point3d A;
	A.input(cin);
	A.output(cout);
	Eigen::MatrixXd m;
	t_mesh::array2matrixd(A, m);
	cout << "tomatrix: \n" << m << endl;
}
void testMesh(igl::opengl::glfw::Viewer &viewer) {
	
	
	//mesh.insert(1.5, 1.5);
	//mesh.insert(2.0, 1.5);
	////mesh.drawControlpolygon(viewer);
	////mesh.drawSurface(viewer, 0.01);
	//mesh.drawTmesh(viewer);
}


int main(int argc,char** argv){
	igl::opengl::glfw::Viewer viewer;
	//testArray();
	//testMesh(viewer);

	//t_mesh::Mesh3d mesh;
	if (mesh.loadMesh("out10.cfg") == 0) {
		cout << "num of nodes: " << mesh.get_num() << endl;
	}
	else {
		cout << "load failing!" << endl;
	}
	//mesh.drawTmesh(viewer);

	
	//mesh.insert(1.5, 1.5);
	//mesh.drawSurface(viewer, 0.01);
	//mesh.drawControlpolygon(viewer);
	
	viewer.callback_key_down = &key_down;
	mesh.draw(viewer, showmesh, showpolygon, showsurface);
	viewer.launch();
	//insert_loop(viewer, mesh);
	//viewer.launch();
    return 0;
}

