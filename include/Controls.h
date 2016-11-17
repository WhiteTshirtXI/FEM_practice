#ifndef _CONTROLS_H_
#define _CONTROLS_H_

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "arcball.h"

#include "Types.h"

namespace Control{

    class ControlBase;

//////////////////////// General Control Functions ///////////

    void initControler(ControlBase* contro);

    void mouseClick(GLFWwindow* window, int button, int action, int mods);

    void mouseMove(GLFWwindow* window, double xpos, double ypos);

    void mouseScroll(GLFWwindow* window, double xoffset, double yoffset);

    void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods);

    void reshape(GLFWwindow* window, int w, int h);

    void help();


//////////////////////// Control Base Virtual Methods ///////////////////////

    /* Control class to handle basic control callback function */
	class ControlBase
	{
	    public:

		    /* controler object initializer initialize */
		    ControlBase(){};

		    /* controler window initialize */
		    virtual void control_init(GLFWwindow* window);

		    /* runtime update */
		    virtual void computeMatrixFromInputs();

		    /* inner variable */
		    virtual glm::vec3 getCamera(){return camera;};
			virtual glm::mat4 getMVP(){ return Projection * View * Model; };
		    virtual glm::mat4 getModel(){return Model;};
		    virtual glm::mat4 getView(){return View;};
		    virtual glm::mat4 getProjection(){return Projection;};

		    /* input control */
		    virtual void mouseClick(GLFWwindow* window, int button, int action, int mods);

		    virtual void mouseMove(GLFWwindow* window, double xpos, double ypos);

		    virtual void mouseScroll(GLFWwindow* window, double xoffset, double yoffset);

		    virtual void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods);

		    virtual void reshape(GLFWwindow* window, int w, int h);

		    virtual void help();

            int shadFlag = 1;

            bool stretchFlag = false;

	    private:
	        /* arcball object */
	        ArcBall arcball;
	        
	        /* rotation quaternion and translation vector for the object */
	        glm::dquat  ObjRot = glm::dquat(1, 0, 0, 0);
			glm::vec3   camera = glm::vec3(0, 0, 2);

	        /* inner variables */
	        int win_width, win_height;
	        double startx, starty;
	        int gButton; 
	        int gState; 

	        /* inner mantained */
	        glm::dmat4 Model;
	        glm::dmat4 View;
	        glm::dmat4 Projection;
	};


}
#endif
