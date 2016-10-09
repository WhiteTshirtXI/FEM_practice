#ifndef _TETRAMESH_H_
#define _TETRAMESH_H_

#include <vector>
#include <string>
#include <Eigen/Dense>

#include "Types.h"

namespace BallonFEM{

    typedef Vertex* pVertex;

    class Vertex
    {
        public:
            Vertex();
            Vertex( Vec3 cord ): m_cord( cord ), m_pos( cord ) {};
            Vertex( Vec3 cord, Vec3 pos ): m_cord( cord ), m_pos( pos ){};

            /* position in material coordinate */
            Vec3 m_cord;
            
            /* position in world space */
            Vec3 m_pos;

            /* index in tetramesh vertex vector */
            int id;

            /* if is fixed */
            bool m_fixed;
    };

    class Face
    {
        public:
            Face(){};
            
            /* compute normal */
            void precomputation();
            
            /* geometry property */
            Vec3 m_normal;

            /* topology property */
            pVertex v[3];
    };

    class Tetra
    {   
        public:
            Tetra(){};

            /* compute volume W and reverse matrix Bm */
            void precomputation();

            /* geometry property */
            Mat3 Bm;    /* reverse matrix of [r1-r4; r2-r4; r3-r4] */
            float W;   /* volume of tetrahedron */

            /* topology property */
            pVertex v[4];
            iVec4 v_id;
    };

    class TetraMesh
    {
        public:
            TetraMesh(){};

            int read( const std::string& filename );
            int write( const std::string& filename );

            /* perform after read to compute surface and tetra properties */
            void precomputation();

            /* vectors of vertices */
            std::vector<Vertex> vertices;

            /* boudary faces of tetramesh */
            std::vector<Face> surface;

            /* tetra mesh */
            std::vector<Tetra> tetrahedrons;

    };

}




#endif //_TETRAMESH_H_
