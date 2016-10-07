#ifndef _CONTROLS_H_
#define _CONTROLS_H_

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

namespace Control{

/* controler initialize */
void control_init(GLFWwindow* window);

/* runtime update */
void computeMatrixFromInputs();

/* inner variable */
glm::vec3 getCamera();
glm::mat4 getMVP();
glm::mat4 getModel();
glm::mat4 getView();
glm::mat4 getProjection();

/* input control */
void mouseClick(GLFWwindow* window, int button, int action, int mods);

void mouseMove(GLFWwindow* window, double xpos, double ypos);

void keyBoard(GLFWwindow* window, int key, int scancode, int action, int mods);

void reshape(GLFWwindow* window, int w, int h);

void help();
}

#endif
