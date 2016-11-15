#include <iostream>
#include <iomanip>

#include "TetraMesh.h"
#include "Dynamic.h"
#include "Control_Simulation.h"
#include "Viewer.h"

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
    BalloonFEM::TetraMesh tetra;
    tetra.read(mesh_name.c_str());

    BalloonFEM::Engine engine(&tetra);
    BalloonFEM::ControlSim controler(&tetra, &engine);

    View::Viewer viewer(&tetra, &controler);

	p_viewer = &(viewer);

    viewer.show();

	return 0;
}
