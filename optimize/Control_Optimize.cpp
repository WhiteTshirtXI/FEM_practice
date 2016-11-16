#include <stdio.h>
#include <ctime>

#include <string>
#include <sstream>
#include <iomanip>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Types.h"
#include "Control_Optimize.h"
#include "Optimizer.h"

namespace BalloonFEM
{

	ControlOpt::ControlOpt(TetraMesh* tetra, Optimizer* engine)
		: ControlBase()
	{
		m_tetra = tetra;
		m_optimizer = engine;
	}

    void ControlOpt::help()
    {  printf("w  -  Wireframe Display\n");
	    printf("f  -  Flat Shading \n");
        printf("r  -  Reset All Parames \n");
        printf("a  -  Add Force \n");
        printf("space - Solve Static Position \n");
	    printf("o  -  Output obj file\n");
	    printf("h  -  Help Information\n");
	    printf("q - quit\n");
    }

    void ControlOpt::keyBoard( GLFWwindow *window, int key, int scancode, int action, int mods)
    {
        if (action != GLFW_PRESS)  return;

		switch (key)
		{
			case GLFW_KEY_F:
				//Flat Shading
				glPolygonMode(GL_FRONT, GL_FILL);
				break;
			case GLFW_KEY_W:
				//Wireframe mode
				glPolygonMode(GL_FRONT, GL_LINE);
				break;
		    case GLFW_KEY_SPACE:
		        Process();
		        break;
			case GLFW_KEY_A:
				AddParameter();
				break;
			case GLFW_KEY_R:
				Reset();
				break;
		    case GLFW_KEY_O:
		        Output();
				break;
			case GLFW_KEY_H:
			case GLFW_KEY_UNKNOWN:
				help();
				break;
			case GLFW_KEY_Q:
			case GLFW_KEY_ESCAPE:
		        glfwSetWindowShouldClose(window, GLFW_TRUE);
				break;
		}

    }

	/* Reset the tetra mesh to original state */
	void ControlOpt::Reset()
	{
		for (size_t i = 0; i < m_tetra->vertices.size(); i++)
		{
			BalloonFEM::Vertex &v = m_tetra->vertices[i];
			v.m_pos = v.m_cord;
			v.m_f_ext = BalloonFEM::Vec3(0);
		}
		force = 0;
		modified = 0;
		shadFlag = 1;
		outputcount = 0;
	}

	void ControlOpt::AddParameter()
	{
	    force += 1e-2;
	    modified = 0;
	    printf("force = %f\n", force);
	}

	/* Do something to tetra mesh */
	void ControlOpt::Process()
	{
		m_optimizer->setCoeff(Vec3(1e4, 10, 1e5));
		m_optimizer->solveOptimal();
		m_optimizer->stepToNext();
		m_optimizer->outputData();
		shadFlag = 1;
	}

	void ControlOpt::Output()
	{
		std::ostringstream info;
		info << "#";
		info << "pressure = ";
		info << m_tetra->holes[0].p << std::endl;

		std::ostringstream ss;
		ss << "output/p_";
		ss << std::setw(7) << std::setfill('0') << outputcount;
		ss << ".obj";
		std::string name(ss.str());

	    printf("Write to %s. \n", name.c_str());

		for (VIter v = m_tetra->vertices.begin(); v != m_tetra->vertices.end(); v++)
			v->m_pos = v->m_cord;
	    m_tetra->write(name, info.str());
		outputcount++;
	}

}

