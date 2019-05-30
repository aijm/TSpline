#include "Test.h"

using namespace std;
using namespace window;


int main(int argc,char** argv){
	/*Test::test_Array();
	Test::test_nlopt();
	Test::test_Integral();*/

	//Test::testBasis();
	//Test::testDerivative();
	
	window::init();
	//Test::test_derOfNurbs();
	//Test::test_lspia();
	Test::testMesh();
	//Test::testSkinning();
	window::launch();
    return 0;
}




