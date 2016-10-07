#ifndef _ARCBALL_H_
#define _ARCBALL_H_

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

namespace {
    
class ArcBall
{
    public:
        ArcBall(){};

        void reset(int width, int height, int xpos, int ypos)
        {
            m_radius = (width > height) ?  width/2 : height/2;

            m_xpos = xpos;
            m_ypos = ypos;

            /* project to unit sphere */
            m_pos = _plane2sphere(xpos, ypos);
        };

        glm::dquat update(int xpos_new, int ypos_new)
        {
            /* project to unit sphere */
            glm::vec3 n_pos = _plane2sphere( xpos_new, ypos_new );
            
            /* actually update 2 times the rotation speed */
            glm::dquat r( glm::dot(m_pos, n_pos) , glm::cross(m_pos, n_pos) );

            m_xpos = xpos_new;
            m_ypos = ypos_new;
            m_pos = n_pos;

            return r;
        };

    private:
        /* project position to unit sphere*/
        glm::vec3 _plane2sphere(int xpos, int ypos)
        {
            double x = xpos / m_radius;
            double y = ypos / m_radius;
            double l =  std::sqrt( x*x + y*y );
            
            /* if curser move out of window */
            if (l>1.0)
                return glm::vec3(x/l, y/l, 0);

            return glm::vec3( x, y, std::sqrt(1 - x*x - y*y) );
        };

        /* radius of arcball */
        double m_radius;

        /* current position of curser */
        double m_xpos, m_ypos;
        

        /* current projection of pos on unit sphere */
        glm::vec3 m_pos;
};

}

#endif
