#include <stdio.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "controls.h"
#include "arcball.h"

#include "Types.h"
#include "Dynamic.h"

extern int shadFlag;

namespace Control{

/* arcball object */
ArcBall arcball;

/* Engine */
BallonFEM::Engine engine;
BallonFEM::TetraMesh* m_tetra;

static double force = 0.01;
static int modified = 0;
void mAddParameter()
{
	force += 0.01;
	modified = 0;
	printf("force = %.2f\n", force);
}

/* Do something to tetra mesh */
void mProcess()
{
	if (modified == 0){
		/* add force to tetra mesh */
		m_tetra->vertices[3].m_f_ext = BallonFEM::Vec3(0.1, 0.1, 0);

		/* specify fixed vertices */
		/*for (size_t i = 0; i < m_tetra->vertices.size(); i++)
		{
			if (abs(m_tetra->vertices[i].m_cord.z) < 1e-2)
				m_tetra->vertices[i].m_fixed = true;
		}*/
		m_tetra->vertices[0].m_fixed = true;
		m_tetra->labelFixedId();

		/* specify rigid bodies */
		std::vector<size_t> rig = {1, 2, 3};
		m_tetra->addRigidBody(rig);

		engine.inputData();
		modified = 1;
	}

    /* compute deformation */
    engine.solveStaticPos();
    engine.stepToNext();
    engine.outputData();

    shadFlag = 1;
}

void mMeasure()
{
	BallonFEM::Vec3 dir(0, 0, force);
	/* add displacement to tetra mesh */
	if (modified == 0){
		for (size_t i = 0; i < m_tetra->vertices.size(); i++)
		{
			BallonFEM::Vertex &v = m_tetra->vertices[i];
			if (abs(v.m_cord.z - 2) < 1e-2)
				v.m_pos = v.m_cord + dir;
		}
		engine.inputData();
		/* specify fixed vertices */
		for (size_t i = 0; i < m_tetra->vertices.size(); i++)
		{
			/* top and bottom faces */
			if (abs(abs(m_tetra->vertices[i].m_cord.z) - 2) < 1e-2)
				m_tetra->vertices[i].m_fixed = true;
		}
		m_tetra->labelFixedId();

		modified = 1;
	}

	/* compute deformation and forces */
	engine.solveStaticPos();
	engine.stepToNext();
	engine.outputData();
	engine.forceTest();

	/* calculate total force */
	/* watch the force on bottom plane is equally*/
	BallonFEM::Vec3 f(0.0);
	for (size_t i = 0; i < m_tetra->vertices.size(); i++)
	{
		BallonFEM::Vertex &v = m_tetra->vertices[i];
		if (abs(v.m_pos.z + 2) < 1e-2)
			f += v.m_velocity;
	}
	printf("Displacement of Z = 2 plane: %f %f %f\n", dir.x, dir.y, dir.z);
	printf("Total Force f = %f, %f, %f \n", f.x, f.y, f.z);
	printf("Young's modulus = %f \n", 4 * f.z / (3.1415926 * dir.z));

	/* re draw */
	shadFlag = 1;
}

/* Reset the tetra mesh to original state */
void mReset()
{
	for (size_t i = 0; i < m_tetra->vertices.size(); i++)
	{
		BallonFEM::Vertex &v = m_tetra->vertices[i];
		v.m_pos = v.m_cord;
		v.m_f_ext = BallonFEM::Vec3(0);
	}
	force = 0;
	modified = 0;
	shadFlag = 1;
}

/* rotation quaternion and translation vector for the object */
glm::dquat  ObjRot(1, 0, 0, 0);
glm::vec3   camera(0, 0, 2);

/* inner variables */
int win_width, win_height;
double startx, starty;
int gButton; 
int gState; 

/* inner mantained */
glm::dmat4 Model;
glm::mat4 View;
glm::mat4 Projection;

/* controler initialize */
void control_init(GLFWwindow* window, BallonFEM::TetraMesh* tetra)
{
    glfwGetWindowSize(window, &win_width, &win_height);
    m_tetra = tetra;
    engine = BallonFEM::Engine(tetra, new BallonFEM::Elastic_neohookean(0.4, 0.4));
}

/* update at each main loop */
void computeMatrixFromInputs()
{
    Model = glm::toMat4(ObjRot);
    View = glm::lookAt(
                        camera,                             /* Camera location */
                        glm::vec3(camera.x, camera.y, 0), /* and looks at this point */
                        glm::vec3(0, 1, 0)                    /* Head is up */
            );
    Projection = glm::perspective(45.0f, win_width/(float) win_height, 0.1f, 100.0f);
}

glm::vec3 getCamera(){
    return camera;
}

glm::mat4 getMVP(){
	return Projection * View * glm::mat4(Model);
}

glm::mat4 getView(){
    return View;
}

glm::mat4 getModel(){
    return glm::mat4(Model);
}

glm::mat4 getProjection(){
    return Projection;
}


/*! mouse click call back function */
void  mouseClick(GLFWwindow* window, int button, int action, int mods) {
	/* set up an arcball around the Eye's center
	switch y coordinates to right handed system  */
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    gState = action;
    
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		gButton = GLFW_MOUSE_BUTTON_LEFT;
		arcball.reset(win_width, win_height, xpos - win_width / 2, win_height - ypos - win_height / 2);
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		startx = xpos;
		starty = ypos;
		gButton = GLFW_MOUSE_BUTTON_RIGHT;
	}
	return;
}

/*! mouse motion call back function */
void mouseMove(GLFWwindow* window, double xpos, double ypos)
{
    glm::vec3   trans;
    glm::dquat  rot;

    if (gState!=GLFW_PRESS)
    {
        return;
    }

	/* rotation, call arcball */
	if (gButton == GLFW_MOUSE_BUTTON_LEFT)
	{
		rot = arcball.update(xpos - win_width / 2, win_height - ypos - win_height / 2);
		ObjRot = rot * ObjRot;
	}

	/*xy translation */
	if (gButton == GLFW_MOUSE_BUTTON_RIGHT)
	{
		double scale = abs( camera.z / win_height );
		trans = glm::vec3(scale*(xpos - startx), scale*(starty - ypos), 0);
		startx = xpos;
		starty = ypos;
		camera = camera - trans;
	}

}

/* mouse middle scroll call back */
void mouseScroll(GLFWwindow* window, double xoffset, double yoffset)
{
	double scale = abs(40. * camera.z / win_height);
	glm::vec3  trans = glm::vec3(0, 0, - scale * yoffset);
	camera = camera + trans;
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
	printf("w  -  Wireframe Display\n");
	printf("f  -  Flat Shading \n");
	printf("s  -  Smooth Shading\n");
	printf("?  -  Help Information\n");
	printf("q - quit\n");
}

/*! Keyboard call back function */
void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS){
        return;
    }

	switch (key)
	{
	case GLFW_KEY_F:
		//Flat Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		break;
	case GLFW_KEY_S:
		//Smooth Shading
		glPolygonMode(GL_FRONT, GL_FILL);
		break;
	case GLFW_KEY_W:
		//Wireframe mode
		glPolygonMode(GL_FRONT, GL_LINE);
		break;
    case GLFW_KEY_SPACE:
        mProcess();
        break;
	case GLFW_KEY_M:
		mMeasure();
		break;
	case GLFW_KEY_A:
		mAddParameter();
		break;
	case GLFW_KEY_R:
		mReset();
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

/*! Called when a "resize" event is received by the window. */
void reshape(GLFWwindow* window, int w, int h)
{
	win_width = w;
	win_height = h;

    /* re project */
    glViewport(0, 0, win_width, win_height);
}

}
