#ifndef _VIEWER_H_
#define _VIEWER_H_

/* this module is used for rendering meshes
 * Mesh type should have at least Vertex and Face element
*/

#include <GLFW/glfw3.h>

#include "TetraMesh.h"
#include "Controls.h"

namespace View{

class Viewer
{
    public:

        /* Mesh type */
        typedef BalloonFEM::TetraMesh Mesh;
		typedef Control::ControlBase Controler;

        Viewer(Mesh* mesh, Controler* controler);

		void refresh(int rebuff = 0);

        void show();

    private:

        /* initialize window and OpenGL context */
        void init_openGL();

        /* setup the drawing */
        void setupGLstate();

        /* prepare program and load shaders */
        void prepareProgram();

        /* buffer data into VBOs*/
        void buffModel();

        /* draw meshes */
        void draw_mesh();

        /* draw axis */
        void draw_axis();

		/* draw force */
		void draw_force();

        /* draw stretches */
        void draw_stretch();

        Mesh* m_mesh;

		Controler* m_controler;

        GLFWwindow* mainWindow;
};

}



#endif
