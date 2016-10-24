#include <iostream>
#include <iomanip>

#include "Viewer.h"
#include "TetraMesh.h"

using namespace std;

View::Viewer *p_viewer;

/*! main function for viewer
*/
int main(int argc, char * argv[])
{
	cout << setprecision(12);

	if (argc < 2)
	{
		cout << "Usage: input.obj" << endl; 
		return -1;
	}

	std::string mesh_name(argv[1]);
		
    /* global mesh */
    BallonFEM::TetraMesh tetra;

    tetra.read(mesh_name.c_str());

    View::Viewer viewer(&tetra);

	p_viewer = &(viewer);

    viewer.show();

	return 0;
}
