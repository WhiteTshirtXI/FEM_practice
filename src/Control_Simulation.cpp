#include <stdio.h>
#include <ctime>

#include <string>
#include <sstream>
#include <iomanip>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Types.h"
#include "Control_Simulation.h"
#include "Dynamic.h"

namespace BalloonFEM
{
    void ControlSim::help()
    {  printf("w  -  Wireframe Display\n");
	    printf("f  -  Flat Shading \n");
        printf("r  -  Reset All Parames \n");
        printf("a  -  Add Force \n");
        printf("space - Solve Static Position \n");
	    printf("o  -  Output obj file\n");
	    printf("h  -  Help Information\n");
	    printf("q - quit\n");
    }

    void ControlSim::keyBoard( GLFWwindow *window, int key, int scancode, int action, int mods)
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
			case GLFW_KEY_M:
		        Measure();
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

	void ControlSim::AddParameter()
	{
	    force += 1e-2;
	    modified = 0;
	    printf("force = %f\n", force);
	}

	/* Do something to tetra mesh */
	void ControlSim::Process()
	{
	    std::vector<double> solveTime;
		//mOutput();
		for (size_t i = 0; i < 1; i++)
		{
			//mAddParameter();
			if (modified == 0){
				/* add disturbution */
				/*for (BalloonFEM::VIter v = m_tetra->vertices.begin(); v != m_tetra->vertices.end(); v++)
				{
					if (v->m_cord.x > 1e-4)
					{
						double theta = 10 * v->m_cord.x;
						v->m_pos = v->m_cord + BalloonFEM::Vec3(0, 0, 0.01*sin(theta));
					}
				}*/

				/* add force to tetra mesh */
				//m_tetra->vertices[3].m_f_ext = BalloonFEM::Vec3(0, -force, 0);
	            /*for (size_t i = 0; i < m_tetra->vertices.size(); i++)
					m_tetra->vertices[i].m_f_ext = BalloonFEM::Vec3(0, -force, 0 );*/

				m_engine->setAirModel(new BalloonFEM::AirModel_Isobaric(force, 0));
				m_engine->inputData();
				modified = 1;
			}

			/* compute deformation */

			std::clock_t start;
			start = std::clock();
			m_engine->solveStaticPos();
			solveTime.push_back((std::clock() - start) / (double) CLOCKS_PER_SEC);

			m_engine->stepToNext();
			m_engine->outputData();
			//BalloonFEM::SpMat K = m_engine->forceTest(force_vec, force_pos);
			//mOutput();
			shadFlag = 1;
		}

		double count = 0;
		for (int i = 0; i < solveTime.size(); i++)
		{
			printf("solve time for step %d: %f s \n", i, solveTime[i]);
			count += solveTime[i];
		}
		printf("total time spent %f s\n", count);
	}

	void ControlSim::Measure()
	{
		AddParameter();
		//m_tetra->vertices[1].m_pos = BalloonFEM::Vec3(cos(force), 0, sin(force));
		for (BalloonFEM::VIter v = m_tetra->vertices.begin(); v != m_tetra->vertices.end(); v++)
		{
			if (v->m_cord.x > 0)
			{
				double theta = force * v->m_pos.x;
				v->m_pos = v->m_cord + BalloonFEM::Vec3((cos(theta) - 1)*v->m_pos.x, sin(theta)*v->m_pos.x, 0);
			}
		}
			
	    m_engine->inputData();
	    //BalloonFEM::SpMat K = m_engine->forceTest(force_vec, force_pos);
		//printf("B location (%.4f, %.4f, %.4f), f[3].z = %.4f\n", cos(force), 0., sin(force), force_vec[3].z);
		shadFlag = 1;
	}

	/* Reset the tetra mesh to original state */
	void ControlSim::Reset()
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

	void ControlSim::Output()
	{
		std::ostringstream info;
		info << "#";
		info << "pressure = ";
		info << force;

		std::ostringstream ss;
		ss << "output/p_";
		ss << std::setw(7) << std::setfill('0') << outputcount;
		ss << ".obj";
		std::string name(ss.str());

	    printf("Write to %s. \n", name.c_str());

	    m_tetra->write(name, info.str());
		outputcount++;
	}

}

