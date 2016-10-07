#include <stdio.h>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "Viewer.h"
#include "controls.h"
using namespace Control;

#include "shaderLoader.h"

#include "MyMesh.h"
using namespace MeshLib;

namespace View{
/* window parameter */
int win_width = 600;
int win_height = 400;

/* GL render */
enum {verVBO, norVBO, faceVBO, cordVBO, num_VBO};
GLuint VAO;
GLuint VBO[num_VBO];
std::vector<glm::vec3> vertex;
std::vector<glm::vec3> norm;
std::vector<unsigned int> faceID;

/* shader data */
GLuint ProgramID;
GLuint MatrixID, CameraID, ModelID, EyeID;

/* render parameter */
glm::vec4 lightColor(0.8, 0.8, 0.8, 1.0);
glm::vec4 globalAmb(0.1, 0.1, 0.1, 1.0);
glm::vec4 eye(1, 0, 0, 0);

/* update eye */
void update_eye()
{
    eye.w += 0.001;
    eye.x = cos(eye.w);
    eye.y = sin(eye.w);
}

/* draw eye */
void draw_eye()
{
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glVertex3d(eye.x, eye.y, eye.z);
    glEnd();
}


Viewer::Viewer(Mesh* mesh)
{
    /* record mesh */
    this->mesh = mesh;

    /* initialize Window and OpenGL */
    init_openGL();

    /* setup drawing environment */
    setupGLstate();

    /* prepare program and shader */
    prepareProgram();

    /* buffer data */
    buff_model();
};

void Viewer::show()
{
    while (!glfwWindowShouldClose(mainWindow))
    {
        /* compute rotation and projection matrix */
        computeMatrixFromInputs();

        /* Render */
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        draw_mesh();
        draw_axis();
        draw_eye();
        
        update_eye();

        /* Swap front and back buffers*/
        glfwSwapBuffers(mainWindow);

        /* Poll for and process event */
        glfwPollEvents();

    }
};

void Viewer::init_openGL()
{
    if (!glfwInit())
    {
        printf("Fail to initialize GLFW.\n");
        return;
    }

    /* create a window and its OpenGL context */
    mainWindow = glfwCreateWindow(win_width, win_height, "Mesh Viewer", NULL, NULL);
    if (!mainWindow)
    {
        printf("Fail to create Window.\n");
        glfwTerminate();
        return;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(mainWindow);

    if (glewInit() != GLEW_OK)
    {
        printf("Fail to initialize GLEW.\n");
        glfwTerminate();
        return;
    }

    /* set interactive callback */
    control_init(mainWindow);
    glfwSetKeyCallback(mainWindow, keyBoard);
    glfwSetMouseButtonCallback(mainWindow, mouseClick);
    glfwSetCursorPosCallback(mainWindow, mouseMove);
    glfwSetWindowSizeCallback(mainWindow, reshape);

};

void Viewer::setupGLstate()
{
	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.35, 0.53, 0.70, 0);
};


void Viewer::prepareProgram()
{
    /* prepare Program */
    ProgramID = shaderLoader("src/vertexShader.glsl", "src/fragmentShader.glsl");
    MatrixID = glGetUniformLocation(ProgramID, "MVP");
    ModelID = glGetUniformLocation(ProgramID, "M");
    CameraID = glGetUniformLocation(ProgramID, "CameraPosition_worldspace");
    EyeID = glGetUniformLocation(ProgramID, "EyePosition_modelspace");
};

void Viewer::buff_model()
{ 
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glGenBuffers(num_VBO, VBO);

    /* prepare vertex & normal data */
    unsigned int N = mesh->numVertices();
    vertex.assign(N, glm::vec3(0,0,0));
    norm.assign(N, glm::vec3(0,0,0));
    unsigned int count = 0;
    for (CMyMesh::_MeshVertexIterator viter(mesh); !viter.end(); ++viter)
    {
        CMyVertex* pV = *viter;
        pV->draw_id() = count;
        CPoint pos = pV->point();
        CPoint nor = pV->normal();
        vertex[count] = glm::vec3(pos[0], pos[1], pos[2]); 
        norm[count] = glm::vec3(nor[0], nor[1], nor[2]);
        count ++;
    }

    glBindBuffer(GL_ARRAY_BUFFER, VBO[verVBO]);
    glBufferData(GL_ARRAY_BUFFER, vertex.size()*sizeof(glm::vec3), &vertex[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, VBO[norVBO]);
    glBufferData(GL_ARRAY_BUFFER, norm.size()*sizeof(glm::vec3), &norm[0], GL_STATIC_DRAW);

    /* prepare faceID data */
    size_t M = mesh->numFaces();
    faceID.assign(M*3, 0);
    count = 0;
    for (CMyMesh::_MeshFaceIterator fiter(mesh); !fiter.end(); ++fiter)
    {
        CMyFace* pF = *fiter;
        CMyHalfEdge* pH = (CMyHalfEdge*)pF->halfedge();
        faceID[3*count] = ((CMyVertex*)pH->source())->draw_id();
        faceID[3*count + 1] = ((CMyVertex*)pH->target())->draw_id();
        faceID[3*count + 2] = ((CMyVertex*)pH->he_next()->target())->draw_id();
        count ++;
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBO[faceVBO]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, faceID.size()*sizeof(unsigned int), &faceID[0], GL_STATIC_DRAW);
};

/*! draw mesh */
void Viewer::draw_mesh()
{
    glUseProgram(ProgramID);
    glBindVertexArray(VAO);

    /* buff uniform variable */
    glm::mat4 MVP = getMVP(); 
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
    glm::mat4 M = getModel();
    glUniformMatrix4fv(ModelID, 1, GL_FALSE, &M[0][0]);
    glm::vec3 camera = getCamera();
    glUniform3f(CameraID, camera.x, camera.y, camera.z);
    glUniform3f(EyeID, eye.x, eye.y, eye.z);
    
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO[verVBO]);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, VBO[norVBO]);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBO[faceVBO]);
    glDrawElements(GL_TRIANGLES, faceID.size(), GL_UNSIGNED_INT, (void*)0);

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
};

void Viewer::draw_axis()
{
    glLineWidth(2.0);
	//x axis
	glColor3f(1.0, 0.0, 0.0);	//red
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(1, 0, 0);
	glEnd();

	//y axis
	glColor3f(0.0, 1.0, 0);		//green
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1, 0);
	glEnd();

	//z axis
	glColor3f(0.0, 0.0, 1.0);	//blue
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 1);
	glEnd();

	glLineWidth(1.0);
};

}
