#include <stdio.h>
#include <stdlib.h>

#include "Viewer.h"
#include "MyMesh.h"

using namespace MeshLib;

/* global mesh */
CMyMesh mesh;

/*! main function for viewer
*/
int main(int argc, char * argv[])
{
	std::cout << std::setprecision(12);

	if (argc < 2)
	{
		std::cout << "Usage: input.m" << std::endl;
		return -1;
	}

	std::string mesh_name(argv[1]);
	if (strutil::endsWith(mesh_name, ".obj"))
	{
		mesh.read_obj(mesh_name.c_str());
	}
	if (strutil::endsWith(mesh_name, ".m"))
	{
		mesh.read_m(mesh_name.c_str());
	}
	if (strutil::endsWith(mesh_name, ".off"))
	{
		mesh.read_off(mesh_name.c_str());
	}

	normalize_mesh(&mesh);
	compute_normal(&mesh);

	mesh.check_gauss_bonnet_theory();

    View::Viewer viewer(&mesh);


    viewer.show();

	return 0;
}
