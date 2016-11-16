#include <iostream>
#include <iomanip>

#include "TetraMesh.h"
#include "Dynamic.h"
#include "Control_Optimize.h"
#include "Viewer.h"

using namespace std;

View::Viewer *p_viewer;

/*! main function for viewer
*/
int main(int argc, char * argv[])
{
	cout << setprecision(12);

	if (argc < 3)
	{
		cout << "Usage: input.vtf target.vtf" << endl; 
		return -1;
	}

	std::string mesh_name(argv[1]);
	std::string target_name(argv[2]);
		
    /* global mesh */
    BalloonFEM::TetraMesh tetra;
    tetra.read(mesh_name.c_str());

	BalloonFEM::TetraMesh target;
	target.read(target_name.c_str());

    BalloonFEM::Optimizer opt(&tetra, &target);
	BalloonFEM::Engine engine(&tetra);
    BalloonFEM::ControlOpt controler(&tetra, &opt, &engine);

    View::Viewer viewer(&tetra, &controler);

	p_viewer = &(viewer);

    viewer.show();

	return 0;
}
