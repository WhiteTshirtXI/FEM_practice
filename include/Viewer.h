#ifndef _VIEWER_H_
#define _VIEWER_H_

/* this module is used for rendering meshes
 * Mesh type should have at least Vertex and Face element
*/

#include <GLFW/glfw3.h>

#include "Tetra.h"

namespace View{

class Viewer
{
    public:

        /* Mesh type */
        typedef Tetra::Surface Mesh;

        Viewer(Mesh* mesh);

        void show();

    private:

        /* initialize window and OpenGL context */
        void init_openGL();

        /* setup the drawing */
        void setupGLstate();

        /* prepare program and load shaders */
        void prepareProgram();

        /* buffer data into VBOs*/
        void buff_model();

        /* draw meshes */
        void draw_mesh();

        /* draw axis */
        void draw_axis();

        Mesh* mesh;

        GLFWwindow* mainWindow;

};

}



#endif
