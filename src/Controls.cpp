#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "Controls.h"
#include "arcball.h"

namespace Control{

////////////////////////// General Control Function ////////////////////////////
    ControlBase* controler;

    void initControler(ControlBase* contro){
        controler = contro;
    }

    void mouseClick(GLFWwindow* window, int button, int action, int mods){
	    controler->mouseClick(window, button, action, mods);
    }

    void mouseMove(GLFWwindow* window, double xpos, double ypos){
	    controler->mouseMove(window, xpos, ypos);
    }

    void mouseScroll(GLFWwindow* window, double xoffset, double yoffset){
	    controler->mouseScroll(window, xoffset, yoffset);
    }

    void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods){
	    controler->keyBoard(window, key, scancode, action, mods);
    }

    void reshape(GLFWwindow* window, int w, int h){
	    controler->reshape(window, w, h);
    }

    void help(){
	    controler->help();
    }

//////////////////////// Control Base Virtual Methods ///////////////////////

    /* controler initialize */
    void ControlBase::control_init(GLFWwindow* mainWindow)
    {
        glfwGetWindowSize(mainWindow, &win_width, &win_height);

    }

    /* update at each main loop */
    void ControlBase::computeMatrixFromInputs()
    {
        Model = glm::toMat4(ObjRot);
        View = glm::lookAt(
                camera,                             /* Camera location */
                glm::vec3(camera.x, camera.y, 0), /* and looks at this point */
                glm::vec3(0, 1, 0)                    /* Head is up */
                );
        Projection = glm::perspective(45.0f, win_width/(float) win_height, 0.1f, 500.0f);
    }


    /*! mouse click call back function */
    void  ControlBase::mouseClick(GLFWwindow* window, int button, int action, int mods) 
    {
	    /* set up an arcball around the Eye's center
	    switch y coordinates to right handed system  */
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        gState = action;
    
	    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	    {
		    gButton = GLFW_MOUSE_BUTTON_LEFT;
		    arcball.reset(win_width, win_height, 
                    xpos - win_width/2, win_height - ypos - win_height/2);
	    }

	    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) 
        {
		    startx = xpos;
		    starty = ypos;
		    gButton = GLFW_MOUSE_BUTTON_RIGHT;
	    }

    }

    /*! mouse motion call back function */
    void ControlBase::mouseMove(GLFWwindow* window, double xpos, double ypos)
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
		    rot = arcball.update(xpos - win_width/2, win_height - ypos - win_height/2);
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
    void ControlBase::mouseScroll(GLFWwindow* window, double xoffset, double yoffset)
    {
	    double scale = abs(40. * camera.z / win_height);
	    glm::vec3  trans = glm::vec3(0, 0, - scale * yoffset);
	    camera = camera + trans;
    }

    /*! helper function to remind the user about commands, hot keys */
    void ControlBase::help()
    {
	    printf("w  -  Wireframe Display\n");
	    printf("f  -  Flat Shading \n");
	    printf("h  -  Help Information\n");
	    printf("q - quit\n");
    }

	/*! Keyboard call back function */
	void ControlBase::keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods)
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
			
			case GLFW_KEY_H:
				help();
				break;
			case GLFW_KEY_Q:
			case GLFW_KEY_ESCAPE:
		        glfwSetWindowShouldClose(window, GLFW_TRUE);
				break;
		}
	}

	/*! Called when a "resize" event is received by the window. */
	void ControlBase::reshape(GLFWwindow* window, int w, int h)
	{
		win_width = w;
		win_height = h;

	    /* re project */
	    glViewport(0, 0, win_width, win_height);
	}

}
