#include <stdio.h>
#include <ctime>

#include <string>
#include <sstream>
#include <iomanip>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Types.h"
#include "Models.h"
#include "Optimizer.h"
#include "Control_Optimize.h"


namespace BalloonFEM
{

	ControlOpt::ControlOpt(TetraMesh* tetra, Optimizer* opt, Engine* engine)
		: ControlBase()
	{
		m_tetra = tetra;
		m_optimizer = opt;
		m_engine = engine;
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
			case GLFW_KEY_S:
				Simulate();
				break;
			case GLFW_KEY_C:
				ChangeAnisoAngle();
				break;
			case GLFW_KEY_M:
				SimulateAniso();
				break;
			case GLFW_KEY_T:
				Target();
				break;
			case GLFW_KEY_V:
                stretchFlag = !stretchFlag;
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
	    force += 0.1;
	    modified = 0;
	    printf("force = %f\n", force);
	}

	/* Do something to tetra mesh */
	void ControlOpt::Process()
	{
		m_optimizer->setCoeff(Vec3(1e1, 1, 1e2));
		m_optimizer->inputData();
		SpVec p(m_tetra->holes.size());
		p << 0.1;
		m_optimizer->setAirPressure(p);

		std::clock_t start;
		start = std::clock();
		m_optimizer->solveOptimalGN();
		printf("Take %.4f solving for optimal\n", (std::clock() - start) / (double)CLOCKS_PER_SEC);

		m_optimizer->stepToNext();
		m_optimizer->outputData();
		shadFlag = 1;
	}

	void ControlOpt::Simulate()
	{
		m_engine->inputData();
		m_engine->solveStaticPos();
		m_engine->stepToNext();
		m_engine->outputData();
		shadFlag = 1;
	}

	void ControlOpt::ChangeAnisoAngle()
	{
		printf("Save current principle stretch direction as material first principle direction \n");
		/* setting material principle direction */
		for (MIter f = m_tetra->films.begin(); f != m_tetra->films.end(); f++)
		for (PIter p = f->pieces.begin(); p != f->pieces.end(); p++)
		{
			p->aniso_angle = p->stretch_angle;
		}
	}

	void ControlOpt::SimulateAniso()
	{
		m_engine->setFilmModel(new Film_aniso_neohookean_3d(0.4, 0.4));
		this->Simulate();
		m_engine->setFilmModel(new Film_neohookean_3d(0.4, 0.4));
	}

	void ControlOpt::Target()
	{
		for (size_t i = 0; i < m_tetra->vertices.size(); i++)
		{
			m_tetra->vertices[i].m_pos = m_optimizer->Target()->vertices[i].m_cord;
		}
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

